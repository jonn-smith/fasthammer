#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <sstream>
#include <utility>
#include <stdexcept>

#include <thread>
#include <future>

#include <assert.h>
#include <x86intrin.h>
#include <string.h>

using namespace std;

#define DO_BATCH_SIMD False

////////////////////////////////////////////////////////////////////////////////

#define LOG(s) do { ::std::cerr << __FILE__ << ":" <<__LINE__ << " : " << (s) << ::std::endl; } while (0);
#define LOGP() do { ::std::cerr << "LOG: POSITION REACHED: " << __FILE__ << ":" <<__LINE__ << ::std::endl; } while (0);
#define LOGSEV(sev, s) do { ::std::cerr << (sev) << ": "; LOG(s); } while (0);
#define LOGMSG(msg, s) do { LOGSEV(msg, s); } while (0);
#define LOGVAR(var) do { ::std::stringstream ss; ss << #var ": "  << var; LOGSEV("DEBUG", ss.str()); } while(0);
#define LOGTOFILE(filename, s) do { ::std::ofstream fout; fout.open( filename ); fout << __FILE__ << ":" << __LINE__ << " : " << (s) << ::std::endl; fout.close(); } while(0);

////////////////////////////////////////////////////////////////////////////////

// From: https://stackoverflow.com/questions/17912933/computing-hamming-distances-to-several-strings-with-sse
static inline int popcnt128(__m128i n) {
  const __m128i n_hi = _mm_unpackhi_epi64(n, n);
  return _mm_popcnt_u64(_mm_cvtsi128_si64(n)) + _mm_popcnt_u64(_mm_cvtsi128_si64(n_hi));
}

// From: https://stackoverflow.com/questions/17912933/computing-hamming-distances-to-several-strings-with-sse
int hamming_vectorized(const unsigned char *p1, unsigned const char *p2, const int len) {
#define MODE (_SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_BIT_MASK | _SIDD_NEGATIVE_POLARITY)
  __m128i smm1 = _mm_loadu_si128 ((__m128i*) p1);
  __m128i smm2 = _mm_loadu_si128 ((__m128i*) p2);
  __m128i ResultMask;

  int iters = len / 16;
  int i;
  int diffs = 0;

  for(i=0; i<iters; i++) {
    ResultMask = _mm_cmpestrm (smm1,16,smm2,16,MODE); 

    diffs += popcnt128(ResultMask);
    p1 = p1+16;
    p2 = p2+16;
    smm1 = _mm_loadu_si128 ((__m128i*)p1);
    smm2 =_mm_loadu_si128 ((__m128i*)p2);
  }

  int mod = len % 16;
  if(mod>0) {
     ResultMask = _mm_cmpestrm (smm1,mod,smm2,mod,MODE); 
     diffs += popcnt128(ResultMask);
  }

  return diffs;
} 

////////////////////////////////////////////////////////////////////////////////

// Adapted From: https://stackoverflow.com/a/33209653
void hammingDistances_SSE(const uint8_t * str, const uint8_t * strings, int * const distances, const int numStrings, const int len)
{
    const int iters = len / 16;

    const __m128i smm1 = _mm_loadu_si128( (__m128i*) str);

    assert((len & 15) == 0);      // m must be a multiple of 16

    for (int j = 0; j < numStrings; j++)
    {
        __m128i smm2 = _mm_loadu_si128( (__m128i*) &strings[ j*(len+1) ]); // len+1, as strings are '\0' terminated

        __m128i diffs = _mm_setzero_si128();

        for (int i = 0; i < iters; i++)
        {
            diffs = _mm_sub_epi8(diffs, _mm_cmpeq_epi8(smm1, smm2));
        }

        diffs = _mm_sad_epu8(diffs, _mm_setzero_si128());
        distances[j] = len - (_mm_extract_epi16(diffs, 0) + _mm_extract_epi16(diffs, 4));
    }
}

vector<size_t> getHamDistBlock_SSE(const vector<string>& whitelist, size_t blockStart, size_t blockEnd) {
	
	vector<size_t> hammingDistances(whitelist[0].length()+1, 0);
	
	size_t d = 0;
	size_t counter = 0;
	
	const size_t numPipsToDisplay = 100;
	const size_t displayThreshold = whitelist.size() / numPipsToDisplay;

	const size_t numJointComparisons = 1024;
	int* jointDistances = new int[numJointComparisons];

	const size_t bufferSize = (whitelist[0].length()+1) * numJointComparisons;
	uint8_t* jointComparisonBuffer = new uint8_t[bufferSize];

	for (size_t i = blockStart; i < whitelist.size() && i < blockEnd; ++i ) {	

		for (size_t j = i+1; j < whitelist.size(); j += numJointComparisons) {
		
			// Get the number of comparisons we're actually performing:
			size_t numActualComparisons = numJointComparisons;
			if ((whitelist.size() - j) < numJointComparisons) {
				numActualComparisons = whitelist.size() - j;
			}

			// Set up our big buffer:
			memset(jointComparisonBuffer, 0, sizeof(char) * bufferSize); 
			memset(jointDistances, 0, sizeof(int) * numJointComparisons); 
			for (size_t k = 0; k < numActualComparisons; ++k) {
				memcpy(jointComparisonBuffer + (sizeof(char) * (k + whitelist[0].length() + 1)), whitelist[j + k].c_str(), whitelist[0].length() + 1);
			}

			// Perform our comparison:
			hammingDistances_SSE( (const uint8_t*) whitelist[i].c_str(), jointComparisonBuffer, jointDistances, numJointComparisons, whitelist[0].length()); 

			// Update our total counts:
			for (size_t k = 0; k < numActualComparisons; ++k) { 
				cerr << jointDistances[k] << " ";
				hammingDistances[jointDistances[k]] = hammingDistances[jointDistances[k]] + 1;
			}
			cerr << endl;
		}

		counter += 1;
		if (counter >= displayThreshold) {
			cerr << "." << flush;
			counter = 0;
		}
	}

	delete[] jointComparisonBuffer;
	delete[] jointDistances;

	for (int i = 0; i < hammingDistances.size(); ++i ){
		cerr << hammingDistances[i] << " ";
	}
	cerr << endl;

	return hammingDistances;
}


////////////////////////////////////////////////////////////////////////////////

size_t hamming(const char* s1, const char* s2, const size_t len) {
	size_t ham = 0;
	for (size_t i = 0 ; i < len; ++i) {
		if (s1[i] != s2[i]) {
			++ham;
		}
	}
	return ham;
}

size_t hamming(const string& s1, const string& s2) {
	size_t ham = 0;
	for (size_t i = 0 ; i < s1.length(); ++i ) {
		if (s1[i] != s2[i]) {
			++ham;
		}
	}
	return ham;
}

void updateMinHamDist(const size_t dist, const size_t index, const vector<string>& whitelist, unordered_map<string, size_t>& minimumHammingDistance) {
	try {
		size_t curDist = minimumHammingDistance.at(whitelist[index]);

		if (dist < curDist) {
			minimumHammingDistance[whitelist[index]] = dist;
		}   
	}   
	catch (const out_of_range& oor) {
		minimumHammingDistance[whitelist[index]] = dist;
	}   
}

pair<vector<size_t>, unordered_map<string, size_t>> getHamDistBlock(const vector<string>& whitelist, size_t blockStart, size_t blockEnd) {
	
	vector<size_t> hammingDistances(whitelist[0].length()+1, 0);
	unordered_map<string, size_t> minimumHammingDistance;

	size_t d = 0;
	size_t counter = 0;
	
	size_t numPipsToDisplay = 100;
	size_t displayThreshold = whitelist.size() / numPipsToDisplay;

	for (size_t i = blockStart; i < whitelist.size() && i < blockEnd; ++i ) {
		for (size_t j = i+1; j < whitelist.size(); ++j) {
			//d = hamming(whitelist[i], whitelist[j]);
			//d = hamming(whitelist[i].c_str(), whitelist[j].c_str(), whitelist[i].length());
			d = hamming_vectorized((const unsigned char*) whitelist[i].c_str(), (const unsigned char*) whitelist[j].c_str(), whitelist[i].length());
			hammingDistances[d] = hammingDistances[d] + 1;
	
			updateMinHamDist(d, i, whitelist, minimumHammingDistance);
			updateMinHamDist(d, j, whitelist, minimumHammingDistance);
		} 

		counter += 1;
		if (counter >= displayThreshold) {
			cerr << "." << flush;
			counter = 0;
		}
	}

	pair<vector<size_t>, unordered_map<string, size_t>> outPair;
	outPair.first = hammingDistances;
	outPair.second = minimumHammingDistance;
	return outPair;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	if (argc < 2) { 
		cerr << argv[0] << " BARCODE_FILE [NUM_THREADS]" << endl;
		return 1;
	} 
	ifstream infile(argv[1]);

	size_t numThreads = thread::hardware_concurrency();
	if (numThreads == 0) {
		numThreads = 2;
	}
	if (argc >= 3) {
		numThreads = stoi(argv[2]); 
	}

	vector<string> whitelist;

	// Read in the barcodes from the whitelist:
	cerr << "Reading in barcodes from " << argv[1] << "..." << endl;
	string bc;
	int len = -1;
	size_t lineNum = 1;
	while (infile >> bc) {
		if (len == -1) {
			len = bc.length();
		}

		if (bc.length() != len) {
			cerr << "ERROR: Encountered barcode of unexpected length on line " << lineNum << "!  Should be of length " << len << ".  Saw: " << bc.length() << " (" << bc << ")" << endl;
			return 1;
		}

		whitelist.push_back(bc);

		lineNum += 1; 
	}
	cerr << "Num barcodes read in: " << whitelist.size() << endl;

	cerr << "Computing hamming distances in " << numThreads << " threads for barcodes of length: " << whitelist[0].length() << endl;

	cerr << "Progress %:" << endl;
	cerr << "        10        20        30        40        50        60        70        80        90         100" << endl;

	vector< future<pair<vector<size_t>, unordered_map<string, size_t>>> > threadOutput;
	const size_t blockSize = whitelist.size() / numThreads;
	for (size_t i = 0; (i*blockSize) < whitelist.size(); ++i) {
#if DO_BATCH_SIMD
		threadOutput.push_back(async(launch::async, getHamDistBlock_SSE, whitelist, i*blockSize, (i+1) * blockSize));
#else
		threadOutput.push_back(async(launch::async, getHamDistBlock, whitelist, i*blockSize, (i+1) * blockSize));
#endif
	}

	// Aggregate outputs:
	vector<size_t> hammingDistances(whitelist[0].length()+1, 0);
	map<string, size_t> minimumHammingDistances;
	for (size_t i = 0; i < threadOutput.size(); ++i) {
		pair<vector<size_t>, unordered_map<string, size_t>> tRet = threadOutput[i].get();

		for (size_t j = 0; j < hammingDistances.size(); ++j) {
			hammingDistances[j] += hammingDistances[j] + tRet.first[j];
		}

		for (auto& v: tRet.second) {
			minimumHammingDistances[v.first] = v.second;
		} 
	}
	
	// Print results:
	cerr << endl << "# Hamming distance distribution:" << endl << endl;
	for (int i = 0; i < hammingDistances.size(); ++i ) {
		cout << i << "\t" << hammingDistances[i] << endl;
	}

	const char* outputFileName = "minimum_hamming_distances.tsv";
	cerr << endl << "Writing per-word min hamming distance to: " << outputFileName << endl;
	ofstream outFile;
	outFile.open(outputFileName);
	if (outFile.is_open()) {
		for (auto& v: minimumHammingDistances) {
			outFile << v.first << "\t" << v.second << endl;
		}
	}
	else {
		cerr << "ERROR: Unable to open output file: " << outputFileName << endl;
	}	
	outFile.close();
	cerr << "DONE" << endl;

	cout << endl;
	return 0;
}

