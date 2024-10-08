// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package util

import "math/bits"

// KmerBaseAt returns the base in pos i (0-based).
func KmerBaseAt(code uint64, k uint8, i uint8) uint8 {
	return uint8(code >> ((k - i - 1) << 1) & 3)
}

// KmerPrefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func KmerPrefix(code uint64, k uint8, n uint8) uint64 {
	return code >> ((k - n) << 1)
}

// KmerSuffix returns the suffix starting from position i (0-based).
// The length of the suffix is k - commonPrefix.
func KmerSuffix(code uint64, k uint8, i uint8) uint64 {
	return code & (1<<((k-i)<<1) - 1)
}

// KmerLongestPrefix returns the length of the longest prefix.
func KmerLongestPrefix(code1, code2 uint64, k1, k2 uint8) uint8 {
	var d uint8
	if k1 >= k2 { // most of the cases
		code1 >>= ((k1 - k2) << 1)
		d = 32 - k2
	} else {
		code2 >>= ((k2 - k1) << 1)
		d = 32 - k1
	}
	return uint8(bits.LeadingZeros64(code1^code2)>>1) - d
}

// MustKmerLongestPrefix returns the length of the longest prefix.
// We assume k1 >= k2.
func MustKmerLongestPrefix(code1, code2 uint64, k1, k2 uint8) uint8 {
	code1 >>= ((k1 - k2) << 1)
	return uint8(bits.LeadingZeros64(code1^code2)>>1) + k2 - 32
}

// KmerHasPrefix checks if a k-mer has a prefix.
func KmerHasPrefix(code uint64, prefix uint64, k1, k2 uint8) bool {
	if k1 < k2 {
		return false
	}
	return code>>((k1-k2)<<1) == prefix
}

// MustKmerHasPrefix checks if a k-mer has a prefix, by assuming k1>=k2.
func MustKmerHasPrefix(code uint64, prefix uint64, k1, k2 uint8) bool {
	return code>>((k1-k2)<<1) == prefix
}

// KmerHasSuffix checks if a k-mer has a suffix.
func KmerHasSuffix(code uint64, suffix uint64, k1, k2 uint8) bool {
	if k1 < k2 {
		return false
	}
	return code&((1<<(k2<<1))-1) == suffix
}

// MustKmerHasSuffix checks if a k-mer has a suffix, by assuming k1>=k2.
func MustKmerHasSuffix(code uint64, suffix uint64, k1, k2 uint8) bool {
	return code&((1<<(k2<<1))-1) == suffix
}

// SharingPrefixKmersMismatch counts the number of mismatch between two k-mers
// sharing with a p-bp prefix.
func SharingPrefixKmersMismatch(code1, code2 uint64, k, p uint8) (n uint8) {
	if p >= k {
		return 0
	}
	var i uint8
	for i = 0; i < k-p; i++ {
		if code1&3 != code2&3 {
			n++
		}
		code1 >>= 2
		code2 >>= 2
	}
	return n
}

// MustSharingPrefixKmersMismatch counts the number of mismatch between two k-mers
// sharing with a p-bp prefix. This function assumes p<k.
func MustSharingPrefixKmersMismatch(code1, code2 uint64, k, p uint8) (n uint8) {
	var i uint8
	for i = 0; i < k-p; i++ {
		if code1&3 != code2&3 {
			n++
		}
		code1 >>= 2
		code2 >>= 2
	}
	return n
}

// SharingPrefixKmersMatches counts the number of matches in the suffix region of two k-mers
// sharing with a p-bp prefix.
func SharingPrefixKmersSuffixMatches(code1, code2 uint64, k, p uint8) (n uint8) {
	if p >= k {
		return 0
	}
	var i uint8
	for i = 0; i < k-p; i++ {
		if code1&3 == code2&3 {
			n++
		}
		code1 >>= 2
		code2 >>= 2
	}
	return n
}

// MustSharingPrefixKmersSuffixMatches counts the number of matches in the suffix region of two k-mers
// sharing with a p-bp prefix.
func MustSharingPrefixKmersSuffixMatches(code1, code2 uint64, k, p uint8) (n uint8) {
	if p >= k {
		return 0
	}
	var i uint8
	for i = 0; i < k-p; i++ {
		if code1&3 == code2&3 {
			n++
		}
		code1 >>= 2
		code2 >>= 2
	}
	return n
}
