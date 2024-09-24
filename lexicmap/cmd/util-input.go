/// Struct for input and output parsing

package cmd

import (
	"strconv"
)

type SearchFields struct {
	query, qseq, sgenome, sseqid, sseq, qcovGnm, hsp, qcovHSP, alenHSP, pident, sstr, cigar, align string
	qlen, gaps, slen, qstart, qend, sstart, send, hits                                             int
}

func SearchFromLine(line string, delimiter byte) SearchFields {
	ncols := 21 // Length of header from other files
	extra_blast_fields := 4
	items := make([]string, ncols)
	stringSplitNByByte(line, delimiter, ncols, &items)
	qlen, _ := strconv.Atoi(items[1])
	hits, _ := strconv.Atoi(items[2])
	gaps, _ := strconv.Atoi(items[10])
	qstart, _ := strconv.Atoi(items[11])
	qend, _ := strconv.Atoi(items[12])
	sstart, _ := strconv.Atoi(items[13])
	send, _ := strconv.Atoi(items[14])
	slen, _ := strconv.Atoi(items[16])

	cigar := ""
	qseq := ""
	sseq := ""
	align := ""

	if len(items) > (ncols - extra_blast_fields) {
		// TODO should work without extra blast fields.
		cigar = items[17]
		qseq = items[18]
		sseq = items[19]
		align = items[20]

	}

	search_output := SearchFields{
		query:   items[0],
		qlen:    qlen,
		hits:    hits,
		sgenome: items[3],
		sseqid:  items[4],
		qcovGnm: items[5],
		hsp:     items[6],
		qcovHSP: items[7],
		alenHSP: items[8],
		pident:  items[9],
		gaps:    gaps,
		qstart:  qstart,
		qend:    qend,
		sstart:  sstart,
		send:    send,
		sstr:    items[15],
		slen:    slen,
		cigar:   cigar,
		qseq:    qseq,
		sseq:    sseq,
		align:   align,
	}

	return search_output
}
