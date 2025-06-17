#include "myutils.h"
#include "syncmerindex.h"
#include "fastaseqsource.h"
#include "seqsource.h"
#include "seqinfo.h"
#include "objmgr.h"



#include <sys/stat.h>
#include <time.h>

off_t get_file_size(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
      return st.st_size;
    }
    return -1;
  }


void cmd_speedbench()
	{
	
	// Defining input file if not given
	const string FileName = {opt(input)} ;
	
	// syncmer parameters
	uint k = optset_k ? opt(k) : 31;
	uint t = optset_t ? opt(t) : 8;
	uint s = k - t; // in basecase it is 31 - 8 = 23
	uint w = 0;  // Not used for closed syncmers, still needed as parameter when building the index

	if (t >= k || t == 0) {
        Die("Invalid parameters: must have 0 < t < k (k=%u, t=%u)", k, t);
    }
	
	// Sequence loading
	FASTASeqSource FSS;
	FSS.Open(FileName);

	SeqInfo *SeqInfo = ObjMgr::GetSeqInfo();
	FSS.GetNext(SeqInfo);
	const byte *Seq = SeqInfo->m_Seq;
	const uint sequence_length = SeqInfo->m_L;



	off_t file_size = get_file_size(FileName.c_str());
    double file_size_mb = file_size / (1024.0 * 1024.0);


	uint type;
	SyncmerType ST = ST_Syncmer2;
	if (optset_algo) {
		type = opt(algo);
    	switch(type) {
			case 1: ST = ST_Syncmer1; break;
			case 2: ST = ST_Syncmer2; break;
			case 3: ST = ST_Syncmer3; break;
			case 4: ST = ST_Syncmer4; break;
			case 5: ST = ST_Syncmer5; break;
			case 6: ST = ST_Syncmer6; break;
			default: ST = ST_Syncmer2; printf("\nValue %u not valid, defaulting to syncmer 2.\n\n", type); type = 2; break;
		}
	}
	
	printf("===================================================\n");
	printf("=== Closed Syncmer (Syncmer %u) Speed Benchmark ===\n", type);
    printf("Parameters: k=%u, t=%u, s=%u\n", k, t, s);
    printf("Input file: %s\n", FileName.c_str());
	printf("Sequence length: %u bp\n", sequence_length);
	printf("File size: %.2f MB\n\n", file_size_mb);
	printf("Syncmer type %u\n", type);

	// SyncmerIndex object creation
	SyncmerIndex Index;

	clock_t start_time = clock();
	Index.Create(ST, k, t, w, Seq, sequence_length);
	clock_t end_time = clock();

	double elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;

	// Getting # k-mers and syncmers for double checking with other implementations
	// uint SyncmerCount = Index.GetSyncmerCount();
    uint KmerCount = Index.GetKmerCount();


	Log("Syncmer computed in %.5lf seconds.\n", elapsed_time) ;
	printf("Results:\n");
    printf("  Time: %.4f seconds\n", elapsed_time);
    printf("  # k-mers: %u\n", KmerCount);
    printf("  # syncmers: %lu (%.2f%%)\n", Index.syncmer_count, 100.0 * Index.syncmer_count / KmerCount);
    printf("  Speed: %.2f MB/s\n", file_size_mb / elapsed_time);

	Index.Clear();

	}