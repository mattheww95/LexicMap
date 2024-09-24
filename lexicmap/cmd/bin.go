// Read binnning module for lexicmap

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
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
const UnspecifiedBin = "NotMapped"

var FastaList []string = []string{".fasta", ".fna", ".fa"}
var FastqList []string = []string{".fastq", ".fq"}

var binCmd = &cobra.Command{
	Use:   "bin",
	Short: "Bin input sequences based on the output of search",
	Long:  "Bin input sequences based on teh output of search",
	Run: func(cmd *cobra.Command, args []string) {
		//opt := getOptions(cmd)
		seq.ValidateSeq = false

		outDirectory := getFlagString(cmd, "out-dir")

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

		_, err = os.Stat(outDirectory)
		if !os.IsNotExist(err) {
			checkError(fmt.Errorf("Output directory should not exist."))
		}

		err = os.Mkdir(outDirectory, 0755)
		checkError(err)

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
				checkError(fmt.Errorf("First line of file is empty.\n"))
			}

			search_val := ProcessInput(scanner.Text())
			searchGenomes[search_val.sgenome] = true
			values_to_append := make([]*SearchFields, 0)
			values_to_append = append(values_to_append, search_val)
			allInputs[search_val.query] = values_to_append
			previous = search_val
		}

		var search_val *SearchFields

		log.Info("Reading binned data.")
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

		// maintain a list of output locations in order to obey the restriction on maximum open files at a time
		outputWrites := make(map[string][]*[]byte)

		log.Info("Assigning queries to genomes.")
		// Organize output files
		for _, file := range files {
			// Initialize outputWrites with each genome to be written too
			for key := range searchGenomes {
				outputWrites[key] = make([]*[]byte, 0, 1)
			}
			outputWrites[UnspecifiedBin] = make([]*[]byte, 0, 1)
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
				if val, ok := allInputs[fastq_id]; ok {
					IdentifyBestHit(val, &outputWrites, record)
					val = nil // remove value from memory
				} else {
					read := record.Format(0)
					outputWrites[UnspecifiedBin] = append(outputWrites[UnspecifiedBin], &read)
					continue
				}
			}
			fastxReader.Close()
			log.Infof("Binning records for %s", file)
			for key, val := range outputWrites {
				var output_file string
				if StringContains(file, &FastqList) {
					output_file = fmt.Sprintf("%s.fastq.gz", key)
				} else if StringContains(file, &FastaList) {
					output_file = fmt.Sprintf("%s.fasta.gz", key)
				} else {
					checkError(fmt.Errorf("Unrecognized input type %s", file))
				}
				output := filepath.Join(outDirectory, output_file)
				outfh, gw, w, err := outStream(output, true, 1)
				checkError(err)
				defer func() {
					outfh.Flush()
					if gw != nil {
						gw.Close()
					}
					w.Close()
				}()

				for _, record := range val {
					outfh.Write(*record)
				}
			}

		}

	},
}

// / Check if any string contains a substring
func StringContains(input string, substr *[]string) bool {
	for _, key := range *substr {
		if strings.Contains(input, key) {
			return true
		}
	}
	return false
}

// / A function for future iteration for identification of an
// / optimal hit per a read if one exists.
func IdentifyBestHit(search_output []*SearchFields, output_genomes *map[string][]*[]byte, fastx *fastx.Record) {
	var previous string = ""
	for _, value := range search_output {
		if value.sgenome != previous {
			read := fastx.Format(0)
			(*output_genomes)[value.sgenome] = append((*output_genomes)[value.sgenome], &read)
		}
		previous = value.sgenome
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

	binCmd.Flags().StringP("out-dir", "o", "binned",
		formatFlagUsage(`Output directory, supports and recommends a ".gz" suffix ("-" for stdout).`))

	binCmd.Flags().StringP("report", "r", "",
		formatFlagUsage(`The generated output of the lexicmap`))

	binCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	binCmd.SetUsageTemplate(usageTemplate(""))
}
