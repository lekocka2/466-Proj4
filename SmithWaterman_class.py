class Smith_Waterman(object):

    """docstring for Smith_Waterman"""
    def __init__(self, seq1,seq2,match=1,mismatch=-1,gap=-1):
        super(Smith_Waterman, self).__init__()
        self.seq1 = seq1
        self.seq2 = seq2
        self.match=match
        self.mismatch=mismatch
        self.gap=gap
        
    def create_score_matrix(self,rows, cols):
        '''Create a matrix of scores representing trial alignments of the two sequences.

        Sequence alignment can be treated as a graph search problem. This function
        creates a graph (2D matrix) of scores, which are based on trial alignments
        of different base pairs. The path with the highest cummulative score is the
        best alignment.
        '''
        score_matrix = [[0 for col in range(cols)] for row in range(rows)]

        # Fill the scoring matrix.
        max_score = 0
        max_pos   = None    # The row and columbn of the highest score in matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                score = self.calc_score(score_matrix, i, j)
                if score > max_score:
                    max_score = score
                    max_pos   = (i, j)

                score_matrix[i][j] = score

        assert max_pos is not None, 'the x, y position with the highest score was not found'

        return score_matrix, max_pos


    def calc_score(self,matrix, x, y):
        '''Calculate score for a given x, y position in the scoring matrix.

        The score is based on the up, left, and upper-left neighbors.
        '''
        similarity = self.match if self.seq1[x - 1] == self.seq2[y - 1] else self.mismatch

        diag_score = matrix[x - 1][y - 1] + similarity
        up_score   = matrix[x - 1][y] + self.gap
        left_score = matrix[x][y - 1] + self.gap

        return max(0, diag_score, up_score, left_score)


    def traceback(self,score_matrix, start_pos):
        '''Find the optimal path through the matrix.

        This function traces a path from the bottom-right to the top-left corner of
        the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
        or both of the sequences being aligned. Moves are determined by the score of
        three adjacent squares: the upper square, the left square, and the diagonal
        upper-left square.

        WHAT EACH MOVE REPRESENTS
            diagonal: match/mismatch
            up:       gap in sequence 1
            left:     gap in sequence 2
        '''

        END, DIAG, UP, LEFT = range(4)
        aligned_seq1 = []
        aligned_seq2 = []
        x, y         = start_pos
        move         = self.next_move(score_matrix, x, y)
        while move != END:
            if move == DIAG:
                aligned_seq1.append(self.seq1[x - 1])
                aligned_seq2.append(self.seq2[y - 1])
                x -= 1
                y -= 1
            elif move == UP:
                aligned_seq1.append(self.seq1[x - 1])
                aligned_seq2.append('-')
                x -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(self.seq2[y - 1])
                y -= 1

            move = self.next_move(score_matrix, x, y)

        aligned_seq1.append(self.seq1[x - 1])
        aligned_seq2.append(self.seq1[y - 1])

        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


    def next_move(self,score_matrix, x, y):
        diag = score_matrix[x - 1][y - 1]
        up   = score_matrix[x - 1][y]
        left = score_matrix[x][y - 1]
        if diag >= up and diag >= left:     # Tie goes to the DIAG move.
            return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
        elif up > diag and up >= left:      # Tie goes to UP move.
            return 2 if up != 0 else 0      # UP move or end.
        elif left > diag and left > up:
            return 3 if left != 0 else 0    # LEFT move or end.
        else:
            # Execution should not reach here.
            raise ValueError('invalid move during traceback')


    def alignment_string(self,aligned_seq1, aligned_seq2):
        '''Construct a special string showing identities, gaps, and mismatches.

        This string is printed between the two aligned sequences and shows the
        identities (|), gaps (-), and mismatches (:). As the string is constructed,
        it also counts number of identities, gaps, and mismatches and returns the
        counts along with the alignment string.

        AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
        ::||::::::||:|::::::: |:  :||:|   <-- alignment string
        CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
        '''
        # Build the string as a list of characters to avoid costly string
        # concatenation.
        idents, gaps, mismatches = 0, 0, 0
        alignment_string = []
        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == base2:
                alignment_string.append('|')
                idents += 1
            elif '-' in (base1, base2):
                alignment_string.append(' ')
                gaps += 1
            else:
                alignment_string.append(':')
                mismatches += 1

        return ''.join(alignment_string), idents, gaps, mismatches


    def print_matrix(self,matrix):
        '''Print the scoring matrix.

        ex:
        0   0   0   0   0   0
        0   2   1   2   1   2
        0   1   1   1   1   1
        0   0   3   2   3   2
        0   2   2   5   4   5
        0   1   4   4   7   6
        '''
        for row in matrix:
            for col in row:
                print('{0:>4}'.format(col))
            print()


    def give_final_result(self):

        # The scoring matrix contains an extra row and column for the gap (-), hence
        # the +1 here.
        rows = len(self.seq1) + 1
        cols = len(self.seq2) + 1

        # Initialize the scoring matrix.
        score_matrix, start_pos = self.create_score_matrix(rows, cols)

        # Traceback. Find the optimal path through the scoring matrix. This path
        # corresponds to the optimal local sequence alignment.
        seq1_aligned, seq2_aligned = self.traceback(score_matrix, start_pos)
        assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'

        # Pretty print the results. The printing follows the format of BLAST results
        # as closely as possible.
        alignment_str, idents, gaps, mismatches = self.alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        # print()
        # print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
                # alength, idents / alength, gaps, alength, gaps / alength))
        # print()
        for i in range(0, alength, 60):
            seq1_slice = seq1_aligned[i:i+60]
            # print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            # print('             {0}'.format(alignment_str[i:i+60]))
            seq2_slice = seq2_aligned[i:i+60]
            # print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            # print()
        
        return seq1_aligned, seq2_aligned, idents / alength, gaps / alength

def main():
    objA = Smith_Waterman("ATGCGCATGCACAAAAGGT","ATGCATGACGGT")
    objB = Smith_Waterman("ATGCAAGCGGTGCATGCACAAAAGGTATGCATGACGGT", "ATGCATGACGGT", 2, -1, -2)
    objA.give_final_result()
    objB.give_final_result()

    
if __name__ == "__main__":
    main()
