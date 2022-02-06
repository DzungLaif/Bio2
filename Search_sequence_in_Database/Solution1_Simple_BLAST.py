# read from file
def read_database(filename):
    f = open(filename)
    db = []
    for line in f:
        db.append(line.rstrip())
    f.close()
    return db

# Hash map to identify all words of size w (a parameter) occurring in query sequence
# The words in the query sequence as key, and a list of occurrence position as values
def build_hash_map(query, width):
    res = {}
    for i in range(len(query)-width+1):
        subseq = query[i:i+width]
    if subseq in res:
        res[subseq].append(i)
    else:
        res[subseq] = [i]
    return res

# find all the matches of words from a sequence with the query
# return a list of hit in tuple type: index of match in query and
# index of match in sequence
def hit_matches (seq, hashedmap, width):
    res = [] # list of tuples
    for i in range(len(seq)-width+1):
        subseq = seq[i:i+width]
        if subseq in hashedmap:
            l = hashedmap[subseq]
            for ind in l:
                res.append((ind, i))
    return res

# from the words in the list, extend the word to both direction
# The result is the tuple: start index in query,start index in sequence
# size of the alignment, and the score
def extends_hit (seq, hit, query, width):
    str_query, str_sequence = hit[0], hit[1]
    #extend foward direction
    matfw = 0
    k = 0
    bestk = 0
    while 2*matfw >= k and str_query+width+k < len(query) and str_sequence+width+k < len(seq):
        if query[str_query+width+k] == seq[str_sequence+width+k]:
            matfw += 1
            bestk = k+1
            k += 1
    size = width + bestk
    # move backwards
    k = 0
    matbw = 0
    bestk = 0
    while 2*matbw >= k and str_query > k and str_sequence > k:
        if query[str_query-k-1] == seq[str_sequence-k-1]:
            matbw += 1
            bestk = k + 1
            k += 1
            size += bestk
    return (str_query-bestk, str_sequence-bestk, size, width+matfw+matbw)

# Find the longest match or best score of each words in the initial list
def hit_best_score(seq, query, hashedmap, width):
    hits = hit_matches(seq, hashedmap, width)
    bestScore = -1.0
    # create a tuple to store results
    best = ()
    for h in hits:
        ext = extends_hit(seq, h, query, width)
        score = ext[3]
        if score > bestScore or (score == bestScore and ext[2] < best[2]):
            bestScore = score
            best = ext
    return best

# Compare all the best results from function above and find the final best result
def best_alignment (db, query, w):
    m = build_hash_map(query, w)
    bestScore = -1.0
    res = (0, 0, 0, 0, 0)
    for k in range(0, len(db)):
        bestSeq = hit_best_score(db[k], query, m, w)
        if bestSeq != ():
            score = bestSeq[3]
            if score > bestScore or (score == bestScore and bestSeq[2] < res[2]):
                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
    if bestScore < 0:
        return ()
    else:
        return res