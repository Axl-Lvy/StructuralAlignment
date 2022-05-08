import numpy as np
import graphviz
from sys import maxsize
import pandas as pd

class Score:
    """Classe d'encapsulation du score et du chemin pour la matrice de Needle
    """
    def __init__(self, score: int, previous: list):
        """Constructeur

        Args:
            score (int): score
            previous (list): liste d'int 1,2 ou 3 selon la direction de la fleche
        """
        self.score = score #the score
        self.previous = previous #the previous pass (only one)
        
    def __str__(self):
        """fonction tostring

        Returns:
            str: le score (pour affichage)
        """
        return str(self.score)
    
    def add(self, a):
        """add a to the previous

        Args:
            a (int): previous, 1 2 ou 3
        """
        self.previous.add(a)

class Needle_Computation:
    """Classe encapsulant toutes les méthodes pour aligner 2 séquences
    """
    
    def __init__(self, seq1: list, seq2: list, id_score=0, sub_score=0, OGP=-4, IGP=-4, blosum=None):
        if isinstance(seq1, str):
            self.seq1 = [seq1]
        else:
            self.seq1 = seq1
        if isinstance(seq2, str):
            self.seq2 = [seq2]
        else:
            self.seq2 = seq2

        if not isinstance(blosum, pd.DataFrame):
            l = [i for i in 'ARNDCQEGHILKMFPSTWYVBZX-']
            blosum = pd.DataFrame([[sub_score]*len(l)]*len(l), columns=l, index=l)
            for i in l:
                blosum[i][i] = id_score
            for i in l:
                blosum[i]['-'] = OGP
                blosum['-'][i] = OGP
            self.blosum = blosum
        else:
            self.blosum = blosum
        self.IGP = IGP
        self.OGP = OGP
        if self.OGP != maxsize:
            self.OGP = min(self.OGP, -self.OGP)
            keys = self.blosum.keys()
            for k in keys:
                self.blosum['-'][k] = self.OGP
                self.blosum[k]['-'] = self.OGP
        self.OGP = self.blosum['-']['A'] #les gap doivent être notés '-' et non '.'
        self.IGP = max(min(self.IGP, -self.IGP), self.OGP)
        X, Y = len(self.seq1[0])+1, len(self.seq2[0])+1
        self.M = np.empty([X, Y], dtype=Score)

    def __str__(self):
        """print the matrix with arrows 
        """
        X,Y = np.shape(self.M)
        s = ""
        for i in range(X):
            for j in range(Y):
                if (2 in self.M[i,j].previous):
                    s+= u"\u2196\t"
                else:
                    s+="\t"
                if (1 in self.M[i,j].previous):
                    s+= u"\u2191\t"
                else:
                    s+="\t"
            s+="\n"
            for j in range(Y):
                if (3 in self.M[i,j].previous):
                    s+=u"\u2190\t"
                else:
                    s+="\t"
                s+=str(self.M[i,j].score)
            s+="\n"
        return s

    def Compute_Matrix(self):
        """Create the needlemanzkefjxh matrix, works with two lists of sequences already aligned

        Returns:
            np.array: needle matrix
        """
        X, Y = len(self.seq1[0])+1, len(self.seq2[0])+1
        self.M[0,0] = Score(0, [])
        for i in range(1,X):
            self.M[i,0] = Score(self.OGP+(i-1)*self.IGP, [1])
        for i in range(1,Y):
            self.M[0,i] = Score(self.OGP+(i-1)*self.IGP, [3])
        for i in range(1, X):
            for j in range(1, Y):
                previous = []
                s1 = self.M[i-1, j].score + (self.IGP if 1 in self.M[i-1,j].previous else self.OGP)
                s3 = self.M[i, j-1].score + (self.IGP if 3 in self.M[i-1,j].previous else self.OGP)
                id_score = 0
                for x in range(len(self.seq1)):
                    for y in range(len(self.seq2)):
                        if self.seq1[x][i-1]=='g' or self.seq2[y][j-1]=='g':
                            print("yoyo")
                        id_score+=self.blosum[self.seq1[x][i-1]][self.seq2[y][j-1]]
                s2 = self.M[i-1, j-1].score + id_score/len(self.seq1)/len(self.seq2)
                s = max([s1, s2, s3])
                if s1 ==s:
                    previous.append(1)
                if s2==s:
                    previous.append(2)
                if s3==s:
                    previous.append(3)
                S = Score(s, previous)
                self.M[i,j]=S
        return self.M
        
    def __add_char_blocks(self, tab: list, char1: list, char2: list)->list:
        """add all the char in char1 and char2 to all athe string in all the alignments in tab
        with tab = ['--AA', 'C-AT', 'CHAT'], ['A-A-', 'C-AT', 'CHAT']] and char1= [-] char2=['T', 'T'] it returns
        ['--AA-', 'C-ATT', 'CHATT'], ['A-A--', 'C-ATT', 'CHATT']]

        Args:
            tab (list): list of alignments
            char1 (list): list of char 1
            char2 (list): list of char 2

        Returns:
            list: new list of sequences
        """
        for i in range(len(char1)):
            tab[i]=char1[i] + tab[i]
        for i in range(len(char2)):
            tab[i+len(char1)]=char2[i] + tab[i+len(char1)]
        return tab

    def Compute_Traceback(self):
        """Compute the traceback, but only one path is considered

        Returns:
            list[str]: the final alignment
        """
        l1, l2 = len(self.seq1), len(self.seq2)
        X, Y = np.shape(self.M)
        X-=1
        Y-=1
        solution = [""]*(l1+l2)
        while(X>0 or Y>0):
            m = self.M[X, Y]
            p = m.previous[0]
            if p==1:
                solution = self.__add_char_blocks(solution, [i[X-1] for i in self.seq1], ["-" for i in range(l2)])
                X-=1
            elif p==2:
                solution = self.__add_char_blocks(solution, [i[X-1] for i in self.seq1], [i[Y-1] for i in self.seq2])
                X-=1
                Y-=1
            elif p==3:
                solution = self.__add_char_blocks(solution, ["-" for i in range(l1)], [i[Y-1] for i in self.seq2])
                Y-=1
        self.alignments = solution
        return self.alignments

class Guide_Tree:


    def __init__(self, seqs: list, blosum: pd.DataFrame, labels=None, IGP=-4, OGP=-4):
        """Constructor

        Args:
            seqs (list): sequences
            blosum (pd.DataFrame): blosum matrix
            labels (list, optional): sequences labels to show on the graph. Defaults to None.
            IGP (int, optional): increase gap penalty. Defaults to -4.
            OGP (int, optional): opening gap penalty. Defaults to -4.
        """

        self.seqs = seqs
        self.blosum = blosum
        if labels==None:
            self.labels = seqs
        else:
            self.labels = labels
        self.IGP = IGP
        self.OGP = OGP
    
    def Compute_Dmatrix(self):
        """Compute the distance matrix

        Returns:
            pd.Dataframe: the distance matrix
        """

        N = len(self.seqs)
        A = np.zeros([N, N])
        for i in range(N):
            for j in range(N):
                if i>j:
                    needle = Needle_Computation(self.seqs[i], self.seqs[j], blosum=self.blosum, IGP=self.IGP, OGP=self.OGP)
                    needle.Compute_Matrix()
                    A[i][j] = needle.M[-1][-1].score
                else:
                    A[i][j] = -maxsize
        for i in range(N):
            for j in range(N):
                if i<j:
                    A[i][j] = A[j][i]
        self.dmatrix = pd.DataFrame(A, index= self.seqs, columns = self.seqs)
        return self.dmatrix

    def __pd_max(df: pd.DataFrame)->tuple:
        """find the max of a pandas dataframe

        Args:
            df (pd.DataFrame): the dataframe

        Returns:
            tuple: (max, i, j) such as df[i][j] = max
        """
        keys = df.keys()
        if len(keys)==1:
            return df[keys[0]][keys[0]], 0, 0
        if len(keys)==0:
            return -maxsize, 0, 0
        M, I, J = -maxsize, 0, 1
        for i in range(len(keys)):
            for j in range(i):
                m = df[keys[i]][keys[j]]
                if m>M:
                    I, J, M = i, j, m
        return M, keys[I], keys[J]

    def Compute_UPGMA(self):
        """Compute the UPGMA method.

        Returns:
            (graphviz.Digraph, list[str]): the graph and the final alignment
        """

        keys = self.labels.copy()
        alignments = {}
        t=0
        if t==0:
            G = graphviz.Digraph()
        for j in range(len(self.dmatrix.keys())):
            G.node(self.dmatrix.keys()[j], label = self.labels[j])
            alignments[self.dmatrix.keys()[j]]=[self.dmatrix.keys()[j]]
            t+=1
        m, i, j = Guide_Tree.__pd_max(self.dmatrix)
        while m<maxsize:
            t+=1
            needle = Needle_Computation(alignments[i], alignments[j], blosum=self.blosum, IGP=self.IGP, OGP=self.OGP)
            needle.Compute_Matrix()
            needle.Compute_Traceback()
            G.node(str(t), label=str(needle.M[-1][-1].score))
            alignments[str(t)] = needle.alignments
            G.edge(i, str(t))
            G.edge(j, str(t))
            l = []
            keys = self.dmatrix.keys()
            for k in keys:
                if k != i or k != j:
                    needle = Needle_Computation(alignments[i], alignments[j], blosum=self.blosum)
                    needle.Compute_Matrix()
                    l.append(needle.M[-1][-1].score)
                else:
                    l.append(-maxsize)
            self.dmatrix.loc[str(t)] = l
            self.dmatrix.insert(len(keys), str(t), l + [maxsize])
            self.dmatrix.drop(i, axis=1, inplace = True)
            self.dmatrix.drop(i, axis=0, inplace = True)
            self.dmatrix.drop(j, axis=1, inplace = True)
            self.dmatrix.drop(j, axis=0, inplace = True)
            m, i, j = Guide_Tree.__pd_max(self.dmatrix)
        self.best_score = needle.M[-1][-1].score
        self.graph = G
        self.consensus = alignments[str(t)]
        return G, alignments[str(t)]
