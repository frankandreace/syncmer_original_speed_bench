#include "myutils.h"
#include "alpha.h"
#include "randseq.h"

static inline byte RandLetter()
	{
	return "ACGT"[randu32()%4];
	}

static inline byte MutateLetter(byte Letter)
	{
	uint r = randu32()%4;
	byte NewLetter = "ACGT"[r];
	if (NewLetter != Letter)
		return NewLetter;
	r = (r + 1)%4;
	NewLetter = "ACGT"[r];
	assert(NewLetter != Letter);
	return NewLetter;
	}

void MakeRandSeq(byte *Seq, uint L)
	{
	for (uint i = 0; i < L; ++i)
		Seq[i] = RandLetter();
	}

void MutateSeq(const byte *InputSeq, uint L, uint PctId, byte *MutatedSeq)
	{
	asserta(PctId <= 100);
	
	uint PctMut = 100 - PctId;
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = InputSeq[i];
		if (randu32()%100 < PctMut)
			Letter = MutateLetter(Letter);
		MutatedSeq[i] = Letter;
		}
	}
