for Sequence in Query do:
    # construct the binary count table
    for Minimiser in Sequence do:
        CountVector = query Minimiser in IBF
        append CountVector to CountTable
    
    for Window in Sequence do:
        # sum over i bitvectors
        # where i = len(sliding window) - k
        WindowCount = sum over slice of CountTable
        
        for BinCount in WindowCount do:
            # compare to k-mer lemma threshold
            if BinCount >= Threshold do:
                add CurrentBin to ResultSet
            
    return ResultSet
