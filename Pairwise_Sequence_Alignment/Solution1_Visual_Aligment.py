import sys

class Draw_Aligment_Plot:
    def __init__(self, seq1, seq2):
        self.sequence1 = seq1
        self.sequence2 = seq2

    # create a matrix that fill with 0
    def create_mat(self, ncolumns, nrows):
        matrix = []
        for i in range(nrows):
            matrix.append([])
            for j in range(ncolumns):
                matrix[i].append(0)
        return matrix

    # create a matrix which the matched character represented by 1 and not match represented by 0
    def dotplot(self):
        matrix_plot = self.create_mat(len(self.sequence1), len(self.sequence2))
        for i in range(len(self.sequence1)):
            for j in range(len(self.sequence2)):
                if self.sequence1[i] == self.sequence2[j]:
                    matrix_plot[i][j] = 1
        return matrix_plot

    # Filter the plot by create a neighborhood window around each position
    # and count the matching characters in this window
    # only fill a neighbor cell if the number of matching character exceeds
    # a limit (stringency)
    def filter_plot(self, window, stringecy):
        matrix_filter = self.create_mat(len(self.sequence1), len(self.sequence2))
        start = int(window/2)
        for i in range(start, len(self.sequence1)-start):
            for j in range(start, len(self.sequence2) - start):
                matches = 0
                l = j - start
                for k in range(i-start, i+start+1):
                    if self.sequence1[k] == self.sequencen2[l]:
                        matches += 1
                    l += 1
                    if matches >= stringecy:
                        matrix_filter[i][j] = 1
        return matrix_filter

    # print plot
    def print_dotplot(self, mat):
        sys.stdout.write(" " + self.sequence2 + "\n")
        for i in range(len(mat)):
            sys.stdout.write(self.sequence1[i])
            for j in range(len(mat[i])):
                if mat[i][j] >= 1:
                    sys.stdout.write("∗")
                else:
                    sys.stdout.write(" ")
            sys.stdout.write("\n")

def test():
    plot_instace = Draw_Aligment_Plot("CGATATACCTAG", "TATGATGGATT")
    matrix_result = plot_instace.filter_plot(5, 4)
    plot_instace.print_dotplot(matrix_result)

test()
