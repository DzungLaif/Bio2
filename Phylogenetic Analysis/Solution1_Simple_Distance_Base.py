from Multiple_Sequence_Alignment import Solution1_Simple_Progressive_MSA_Part_I

# Create a binary tree with the input sequences in the leaves
# and the internal nodes is bifurcation points - mutation events
# internal node contain ít height, leaves contain the index of the sequences
class BinaryTree:
    def __init__(self, val, dist=0, left=None, right=None):
        # integer - an index of the sequence in the leaves, -1 in internal node
        self.value = val
        # height of the node, (0 at the leaves)
        self.distance = dist
        # left sub tree
        self.left = left
        # left right tree
        self.right = right

    def print_tree(self):
        self.print_tree_rec(0, "Root")

    def print_tree_rec(self, level, side):
        tabs = ""
        for i in range(level):
            tabs += "\t"
        if self.value >= 0:
            print(tabs, side, "− value:", self.value)
        else:
            print(tabs, side, "− Dist.: ", self.distance)
            if (self.left != None):
                self.left.print_tree_rec(level + 1, "Left")
            if (self.right != None):
                self.right.print_tree_rec(level + 1, "Right")

    def get_cluster(self):
        res = []
        if self.value >= 0:
            res.append(self.value)
        else:
            if (self.left != None):
                res.extend(self.left.get_cluster())
            if (self.left != None):
                res.extend(self.right.get_cluster())
        return res

# UPGMA Algorithm
# Step 1: create a nummatrix (distance matrix) to store and manipulate.
# this will be the input matrix for UPGMA class
class NumMatrix:
    def __init__(self, rows, cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for j in range(cols):
                self.mat[i].append(0.0)

    def __getitem__(self, n):
        return self.mat[n]

    def num_rows(self):
        return len(self.mat)

    def num_cols(self):
        return len(self.mat[0])

    def get_value(self, i, j):
        if i > j: return self.mat[i][j]
        else: return self.mat[j][i]

    def set_value(self, i, j, value):
        if i > j: self.mat[i][j] = value
        else: self.mat[j][i] = value

    def print_mat(self):
        for r in self.mat: print(r)
        print()

    def min_dist_indexes(self):
        m = self.mat[1][0]
        res = (1, 0)
        for i in range(1, self.num_rows()):
            for j in range(i):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
                    res = (i, j)
        return res

    def add_row(self, newrow):
        self.mat.append(newrow)

    def add_col(self, newcol):
        for r in range(self.num_rows()):
            self.mat[r].append(newcol[r])

    def remove_row(self, ind):
        del self.mat[ind]

    def remove_col(self, ind):
        for r in range(self.num_rows()):
            del self.mat[r][ind]

    def copy(self):
        newm = NumMatrix(self.num_rows(), self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                newm.mat[i][j] = self.mat[i][j]
        return newm

# Step 3: Create Hierarchical Clustering
class HierarchicalClustering:
    def __init__(self, matx_of_distance: NumMatrix):
        self.matdists = matx_of_distance

    def execute_clustering(self):
        # initialization of the tree leaves and matrix
        trees = []
        for i in range(self.matdists.num_rows()):
            t = BinaryTree(i)
            trees.append(t)
        tableDist = self.matdists.copy()
        for k in range(self.matdists.num_rows(), 1, -1):
            ## minimum distance in D
            mins = tableDist.min_dist_indexes()
            i, j = mins[0], mins[1]
            ## create new tree joining clusters
            n = BinaryTree(-1, tableDist.get_value(i, j) / 2.0, trees[i], trees[j])
            if k > 2:
                ## remove trees being joined from the list
                ti = trees.pop(i)
                tj = trees.pop(j)
                ## calculating distances for new cluster
                dists = []
                for x in range(tableDist.num_rows()):
                    if x != i and x != j:
                        si = len(ti.get_cluster())
                        sj = len(tj.get_cluster())
                        d = (si*tableDist.get_value(i, x) + sj*tableDist.get_value(j, x)) / (si + sj)
                        dists.append(d)
                ## updating the matrix
                tableDist.remove_row(i)
                tableDist.remove_row(j)
                tableDist.remove_col(i)
                tableDist.remove_col(j)
                tableDist.add_row(dists)
                tableDist.add_col([0]*(len(dists) + 1))
                ## add the new tree to the set to handle
                trees.append(n)
            else: return n

# Step 3: Implement UPGMA
class UPGMA:
    def __init__(self, seqs: [], alseq: Solution1_Simple_Progressive_MSA_Part_I.PairwiseAlignment):
        self .seqs = seqs
        self .alseq = alseq
        self .create_mat_dist()

    def create_mat_dist(self):
        self.matdist = NumMatrix(len(self.seqs), len(self.seqs))
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needleman_Wunsh_Algorithsm(s1, s2)
                alin = self.alseq.recover_align()
                ncd = 0
                for k in range(len(alin)):
                    col = alin.column(k)
                    if col[0] != col[1]: ncd += 1
                self.matdist.set_value(i, j, ncd)

    def run(self):
        ch = HierarchicalClustering(self .matdist)
        t = ch.execute_clustering()
        return t

def test():
    seq1 = Solution1_Simple_Progressive_MSA_Part_I.MySeq("ATAGCGAT")
    seq2 = Solution1_Simple_Progressive_MSA_Part_I.MySeq("ATAGGCCT")
    seq3 = Solution1_Simple_Progressive_MSA_Part_I.MySeq("CTAGGCCC")
    seq4 = Solution1_Simple_Progressive_MSA_Part_I.MySeq("CTAGGCCT")
    sm = Solution1_Simple_Progressive_MSA_Part_I.SubstMatrix()
    sm.create_submat(1, -1, "ACGT")
    alseq = Solution1_Simple_Progressive_MSA_Part_I.PairwiseAlignment(sm, -2)
    up = UPGMA([seq1, seq2, seq3, seq4], alseq)
    arv = up.run()
    arv.print_tree()