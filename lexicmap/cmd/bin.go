// Read binnning module for lexicmap

package cmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"
	"time"
	"unsafe"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

const ShortHeader string = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen"
const LongHeader string = "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tsstr\tslen\tcigar\tqseq\tsseq\talign"
const UnspecifiedBin = "NotMapped"
const UniqueBinned string = "Unique"
const AllBinned string = "All"

var FastaList []string = []string{".fasta", ".fna", ".fa"}
var FastqList []string = []string{".fastq", ".fq"}

var binCmd = &cobra.Command{
	Use:   "bin",
	Short: "Bin input sequences based on the output of search",
	Long:  "Bin input sequences based on teh output of search",
	Run: func(cmd *cobra.Command, args []string) {
		opt := getOptions(cmd)
		seq.ValidateSeq = false

		var fhLog *os.File
		if opt.Log2File {
			fhLog = addLog(opt.LogFile, opt.Verbose)
		}

		verbose := opt.Verbose
		outputLog := opt.Verbose || opt.Log2File

		timeStart := time.Now()
		defer func() {
			if outputLog {
				log.Info()
				log.Infof("elapsed time: %s", time.Since(timeStart))
				log.Info()
			}
			if opt.Log2File {
				fhLog.Close()
			}
		}()

		bin_unique := getFlagBool(cmd, "bin-unique-reads")
		outDirectory := getFlagString(cmd, "out-dir")

		bufferSizeS := getFlagString(cmd, "buffer-size")
		if bufferSizeS == "" {
			checkError(fmt.Errorf("value of buffer size. supported unit: K, M, G"))
		}

		files := getFileListFromArgsAndFile(cmd, args, true, "infile-list", true)
		if len(files) == 1 {
			if isStdin(files[0]) {
				checkError(fmt.Errorf("  no files given, cant read from stdin"))
			} else {
				log.Infof("  %d input file given: %s", len(files), files[0])
			}
		} else {
			log.Infof("  %d input file(s) given", len(files))
		}

		var err error
		bufferSize, err := ParseByteSize(bufferSizeS)
		if err != nil {
			checkError(fmt.Errorf("invalid value of buffer size. supported unit: K, M, G"))
		}

		//outputLog := opt.Verbose || opt.Log2File
		report := getFlagString(cmd, "report")

		if isStdin(report) {
			log.Info("  no files given, reading from stdin")
		} else if outputLog {
			log.Infof("Input file given: %s", report)
		}

		_, err = os.Stat(outDirectory)
		if !os.IsNotExist(err) {
			checkError(fmt.Errorf("output directory should not exist"))
		}

		err = os.Mkdir(outDirectory, 0755)
		log.Infof("Created output directory: %s", outDirectory)
		checkError(err)

		buf := make([]byte, bufferSize)
		var fh *xopen.Reader
		var line string
		var scanner *bufio.Scanner
		var blastOutput bool = false
		log_read_interval := 1000

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
				checkError(fmt.Errorf("first line of file is empty"))
			}

			search_val := ProcessInput(scanner.Text())
			searchGenomes[search_val.sgenome] = true
			values_to_append := make([]*SearchFields, 0)
			values_to_append = append(values_to_append, search_val)
			allInputs[search_val.query] = values_to_append
			previous = search_val
		}

		var search_val *SearchFields

		if outputLog {
			log.Info("Reading binned data.")
		}

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
		allInputs_l := len(allInputs)
		if outputLog {
			log.Infof("Sequences listed in output: %d", allInputs_l)
		}

		checkError(scanner.Err())
		checkError(fh.Close())
		var record *fastx.Record
		SequenceTracker := 0
		var BasesInMemory uint64 = 0
		var FlushBasesThreshold uint64 = 4_000_000_000 // Flush when x amount of bases are in memory

		log.Info("Assigning queries to genomes.")
		// Organize output files
		for _, file := range files {
			log.Infof("Processing: %s", file)

			// Initialize outputWrites with each genome to be written too
			outputWrites := CreateOutputGroups(&searchGenomes, UnspecifiedBin)
			outputWrites_unq := CreateOutputGroups(&searchGenomes, UnspecifiedBin)
			fastxReader, err := fastx.NewReader(nil, file, "")

			checkError(err)
			timeStart1 := time.Now()
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
				BasesInMemory = BasesInMemory + uint64(len(record.Seq.Seq))
				read := record.Format(0)
				if val, ok := allInputs[fastq_id]; ok {
					IdentifyBestHit(val, &outputWrites, &outputWrites_unq, &read, bin_unique)
					val = nil // remove value from memory
				} else {
					// Unspecified is always unique, but no need to track it twice
					outputWrites[UnspecifiedBin] = Append(outputWrites[UnspecifiedBin], &read)
				}
				SequenceTracker++
				if (SequenceTracker%log_read_interval) == 0 && verbose && outputLog {
					speed := float64(SequenceTracker) / time.Since(timeStart1).Minutes()
					fmt.Fprintf(os.Stderr, "Processed: %d of %d records %.3f matches per minute \r", SequenceTracker, allInputs_l, speed)
					if BasesInMemory >= FlushBasesThreshold {
						// Flush the outputs periodically to prevent the maps from growing too large
						// causing paging to disk
						if bin_unique {
							WriteBinnedReads(&outputWrites_unq, file, outDirectory, UniqueBinned, false, opt.CompressionLevel)
							WriteBinnedReads(&outputWrites, file, outDirectory, AllBinned, true, opt.CompressionLevel)
						} else {
							WriteBinnedReads(&outputWrites, file, outDirectory, "", true, opt.CompressionLevel)
						}
						BasesInMemory = 0
					}
				}
			}

			fastxReader.Close()
			if outputLog {
				log.Infof("Binning records for %s", file)
			}

			if bin_unique {
				WriteBinnedReads(&outputWrites_unq, file, outDirectory, UniqueBinned, false, opt.CompressionLevel)
				WriteBinnedReads(&outputWrites, file, outDirectory, AllBinned, true, opt.CompressionLevel)
			} else {
				WriteBinnedReads(&outputWrites, file, outDirectory, "", true, opt.CompressionLevel)
			}
			BasesInMemory = 0
		}

		if outputLog {
			log.Info("Done")
		}

	},
}

// / Needed to prevent the overhead of string creation as map keys need to be strings
// / in golang. This function comes from: https://syslog.ravelin.com/byte-vs-string-in-go-d645b67ca7ff
func ByteSliceToString(bs []byte) string {
	return *(*string)(unsafe.Pointer(&bs))
}

func CreateOutputGroups(genomes *map[string]bool, unspecified_bin string) map[string][]*[]byte {
	outputWrites := make(map[string][]*[]byte)
	for key := range *genomes {
		outputWrites[key] = make([]*[]byte, 0, 10)
	}
	outputWrites[unspecified_bin] = make([]*[]byte, 0, 10)
	return outputWrites
}

// / Get output file name per a file
func GetOutputFile(input_file string, output_directory string, output_name string, nested_directory string) string {
	var output_file string
	var output string
	if StringContains(input_file, &FastqList) {
		output_file = fmt.Sprintf("%s.fastq.gz", output_name)
	} else if StringContains(input_file, &FastaList) {
		output_file = fmt.Sprintf("%s.fasta.gz", output_name)
	} else {
		checkError(fmt.Errorf("unrecognized input type %s", input_file))
	}
	if nested_directory == "" {
		output = filepath.Join(output_directory, output_file)
	} else {
		output = filepath.Join(output_directory, nested_directory, output_file)
	}
	return output
}

func WriteBinnedReads(outputs *map[string][]*[]byte, file_name string, output_directory string, nested_string string, clear_records bool, compression_level int) {
	for key, val := range *outputs {
		output := GetOutputFile(file_name, output_directory, key, nested_string)
		outfh, gw, w, err := outStream(output, strings.HasSuffix(output, ".gz"), compression_level)
		checkError(err)

		for _, record := range val {
			outfh.Write(*record)
			if clear_records {
				record = nil
			}
		}
		if clear_records {
			val = nil
		}
		outfh.Flush()
		if gw != nil {
			gw.Close()
		}
		w.Close()
	}
	runtime.GC()
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
func IdentifyBestHit(search_output []*SearchFields, output_genomes *map[string][]*[]byte, output_genomes_unq *map[string][]*[]byte, read *[]byte, id_best bool) {
	var previous string = ""
	sgenomes_out := 0
	for _, value := range search_output {
		if value.sgenome != previous {
			sgenomes_out++
			(*output_genomes)[value.sgenome] = Append((*output_genomes)[value.sgenome], read)
		}
		previous = value.sgenome
	}

	if sgenomes_out == 1 && id_best {
		(*output_genomes_unq)[previous] = Append((*output_genomes_unq)[previous], read)
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
	return match
}

func Append(slice []*[]byte, new_value ...*[]byte) []*[]byte {
	n := len(slice)
	total := len(slice) + len(new_value)
	if n == cap(slice) {
		newSize := total * 2 // grow array by 2
		newSlice := make([]*[]byte, total, newSize)
		copy(newSlice, slice)
		slice = newSlice
	}
	slice = slice[:total]
	copy(slice[n:], new_value)
	return slice
}

func init() {
	utilsCmd.AddCommand(binCmd)

	binCmd.Flags().StringP("out-dir", "o", "binned",
		formatFlagUsage(`Output directory, supports and recommends a ".gz" suffix ("-" for stdout).`))

	binCmd.Flags().StringP("report", "r", "",
		formatFlagUsage(`The generated output of the lexicmap`))

	binCmd.Flags().StringP("buffer-size", "b", "20M",
		formatFlagUsage(`Size of buffer, supported unit: K, M, G. You need increase the value when "bufio.Scanner: token too long" error reported`))

	binCmd.Flags().BoolP("bin-unique-reads", "u", true,
		formatFlagUsage("Create separate reads from unique source into a separate folder."))

	binCmd.SetUsageTemplate(usageTemplate(""))
}
