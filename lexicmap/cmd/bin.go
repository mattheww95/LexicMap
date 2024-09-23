// Read binnning module for lexicmap

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"path/filepath"
	"regexp"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

const ShortHeader string = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen"
const LongHeader string = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tcigar\tqseq\tsseq\talign"

var binCmd = &cobra.Command{
	Use:   "bin",
	Short: "Bin input sequences based on the output of search",
	Long:  "Bin input sequences based on teh output of search",
	Run: func(cmd *cobra.Command, args []string) {
		//opt := getOptions(cmd)
		seq.ValidateSeq = false

		outFile := getFlagString(cmd, "out-file")

		bufferSizeS := getFlagString(cmd, "buffer-size")
		if bufferSizeS == "" {
			checkError(fmt.Errorf("value of buffer size. supported unit: K, M, G"))
		}
		var err error
		bufferSize, err := ParseByteSize(bufferSizeS)
		if err != nil {
			checkError(fmt.Errorf("invalid value of buffer size. supported unit: K, M, G"))
		}

		//verbose := opt.Verbose
		//outputLog := opt.Verbose || opt.Log2File
		report := getFlagString(cmd, "report")

		if isStdin(report) {
			log.Info("  no files given, reading from stdin")
		} else {
			log.Infof("Input file given: %s", report)
		}

		outFileClean := filepath.Clean(outFile)
		if !isStdin(report) && filepath.Clean(report) == outFileClean {
			checkError(fmt.Errorf("out file should not be one of the input file"))
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if len(files) == 1 {
			if isStdin(files[0]) {
				log.Info("  no files given, reading from stdin")
			} else {
				log.Infof("  %d input file given: %s", len(files), files[0])
			}
		} else {
			log.Infof("  %d input file(s) given", len(files))
		}

		buf := make([]byte, bufferSize)
		var fh *xopen.Reader
		var line string
		var scanner *bufio.Scanner
		var blastOutput bool = false

		//var value string
		fh, err = xopen.Ropen(report)
		checkError(err)

		scanner = bufio.NewScanner(fh)
		scanner.Buffer(buf, int(bufferSize))
		scanner.Scan()
		headerLine := false

		// Check if additional blast fields are used
		blastOutput = CheckRegex(scanner.Text(), LongHeader)
		if CheckRegex(scanner.Text(), ShortHeader) || blastOutput {
			headerLine = true
		}

		allInputs := make(map[string][]*SearchFields)
		searchGenomes := make(map[string]bool) // Making a set of target genomes
		var previous *SearchFields
		if headerLine {
			scanner.Scan()
			if scanner.Text() == "" {
				checkError()
			}

			search_val := ProcessInput(scanner.Text())
			searchGenomes[search_val.query] = true
			values_to_append := make([]*SearchFields, 0)
			values_to_append = append(values_to_append, search_val)
			allInputs[search_val.query] = values_to_append
			previous = search_val
		}

		var search_val *SearchFields

		for scanner.Scan() {
			line = scanner.Text()
			if line == "" {
				continue
			}
			search_val = ProcessInput(line)
			if _, ok := searchGenomes[search_val.sgenome]; !ok {
				searchGenomes[search_val.sgenome] = true
			}

			if previous.query != search_val.query { // Add new value if read is new or previous is nil
				values_to_append := make([]*SearchFields, 0)
				values_to_append = append(values_to_append, search_val)
				allInputs[search_val.query] = values_to_append
				previous = search_val
			} else {
				allInputs[search_val.query] = append(allInputs[search_val.query], search_val)
			}
		}

		checkError(scanner.Err())
		checkError(fh.Close())
		var record *fastx.Record
		for _, file := range files {
			fastxReader, err := fastx.NewReader(nil, file, "")
			checkError(err)
			for {
				record, err = fastxReader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					checkError(err)
					break
				}
				fastq_id := string(record.ID)
				//if val, ok := allInputs[fastq_id]; ok {
				if _, ok := allInputs[fastq_id]; ok {
					//IdentifyBestHit(val)
					continue
				} else {
					//log.Infof("Missing %s in search output.", fastq_id)
					break
				}

			}
			fastxReader.Close()
		}
		fmt.Println(searchGenomes)
	},
}

// / A function for future iteration for identification of an
// / optimal hit per a read if one exists.
func IdentifyBestHit(search_output []*SearchFields) {
	for _, value := range search_output {
		fmt.Printf("%+v\n", value)
	}
}

// / Process new input line
func ProcessInput(line string) *SearchFields {
	line_stripped := strings.TrimRight(line, "\r\n")
	new_line := SearchFromLine(line_stripped, '\t')
	return &new_line
}

// / Determine if the header field is present in the file input
// / line string: The incoming text line to check
// / header_match string: The pattern to test for
func CheckRegex(line string, header_match string) bool {
	match, err := regexp.MatchString(header_match, line)
	checkError(err)

	if match {
		return true
	}
	return false
}

func init() {
	utilsCmd.AddCommand(binCmd)

	binCmd.Flags().StringP("out-file", "o", "-",
		formatFlagUsage(`Out file, supports and recommends a ".gz" suffix ("-" for stdout).`))

	binCmd.Flags().StringP("report", "r", "",
		formatFlagUsage(`The generated output of the lexicmap`))

	binCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	binCmd.SetUsageTemplate(usageTemplate(""))
}
