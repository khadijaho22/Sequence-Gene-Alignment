from struct import pack
import numpy as np
from tkinter import *
from Bio import Align
from tkinter import filedialog
from tkinter import messagebox
from PIL import ImageTk, Image

list=[]


def smith_waterman(seq1, seq2, match_score=2, mismatch_score=-1, gap_penalty=-2):
 
   
    # Initialize the scoring matrix with zeros
    nrow = len(seq1) + 1
    ncol = len(seq2) + 1
    score_matrix = np.zeros((nrow, ncol))

    # Initialize the highest score and index
    max_score = 0
    max_i, max_j = 0, 0

    # Fill in the scoring matrix by considering three possible ways to get to each cell
    for i in range(1, nrow):
        for j in range(1, ncol):
            match = score_matrix[i-1,j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            gap1 = score_matrix[i-1,j] + gap_penalty
            gap2 = score_matrix[i-1,j-1] + gap_penalty
            score_matrix[i,j] = max(0, match, gap1, gap2)

            # Update the highest score and index
            if score_matrix[i,j] > max_score:
                max_score = score_matrix[i,j]
                max_i, max_j = i, j


    align1 = ''
    align2 = ''
    i, j = max_i, max_j
    while score_matrix[i,j] != 0:
        score = score_matrix[i,j]
        score_diag = score_matrix[i-1,j-1]
        score_up = score_matrix[i-1,j]
        score_left = score_matrix[i,j-1]
        if score == score_diag + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score == score_up + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif score == score_left + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, max_score

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):

    # Initialize the scoring matrix with zeros
    nrow = len(seq1) + 1
    ncol = len(seq2) + 1
    score_matrix = np.zeros((nrow, ncol))

    # Initialize the first row and column with gap penalty
    for i in range(1, nrow):
        score_matrix[i,0] = i * gap_penalty
    for j in range(1, ncol):
        score_matrix[0,j] = j * gap_penalty

    # Fill in the scoring matrix by considering three possible ways to get to each cell
    for i in range(1, nrow):
        for j in range(1, ncol):
            match = score_matrix[i-1,j-1] + (match_score if seq1[i-1] ==seq2[j-1] else mismatch_score)
            gap1 = score_matrix[i-1,j] + gap_penalty
            gap2 = score_matrix[i,j-1] + gap_penalty
            score_matrix[i,j] = max(match, gap1, gap2)

    # Backtrack to find the optimal alignment
    align1 = ''
    align2 = ''
    i = nrow - 1
    j = ncol - 1
    while i > 0 and j > 0:
        score = score_matrix[i,j]
        score_diag = score_matrix[i-1,j-1]
        score_up = score_matrix[i-1,j]
        score_left = score_matrix[i,j-1]
        if score == score_diag + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score == score_up + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif score == score_left + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = '-' + align2
        i -= 1
    while j > 0:
        align1 = '-' + align1
        align2 = seq2[j-1] + align2
        j -= 1

    return align1, align2, score_matrix[nrow-1, ncol-1]

def global_alignment_brute(seq1, seq2,match_score=1, mismatch_score=-1, gap_penalty=-1):
    # Generate all possible alignments recursively
     
    all_alignments = []
    def generate_alignments(align1='', align2='', i=0, j=0):
        if i == len(seq1) and j == len(seq2):
            all_alignments.append((align1, align2))
        if i < len(seq1):
            generate_alignments(align1 + seq1[i], align2 + '-', i+1, j)
        if j < len(seq2):
            generate_alignments(align1 + '-', align2 + seq2[j], i, j+1)
        if i < len(seq1) and j < len(seq2):
            generate_alignments(align1 + seq1[i], align2 + seq2[j], i+1, j+1)
    
    generate_alignments()
    
    # Compute the score for each alignment
    aligned_scores = []
    for align1, align2 in all_alignments:
        score = 0
        for c1, c2 in zip(align1, align2):
            if c1 == '-' or c2 == '-':
                score += gap_penalty
            elif c1 == c2:
                score += match_score
            else:
                score += mismatch_score
        aligned_scores.append(score)
    
    # Find the alignment with the highest score
    max_score_index = aligned_scores.index(max(aligned_scores))
    align1, align2 = all_alignments[max_score_index]
    
    # print alignment and its score
    print(align1, align2, max(aligned_scores))
    
def greedy_align(seq11, seq2):
    aligned = ""
    seq1 = seq11
    len_diff = abs(len(seq1) - len(seq2))
    max_len = max(len(seq1), len(seq2))
    if len(seq1) > len(seq2):
      seq1 = seq1[0:(len(seq1)-len_diff)]
    else:
       seq2 = seq2[0:(len(seq2)-len_diff)]
    for i in range(len(seq1)):
      if seq1[i] == seq2[i]:
        aligned += seq1[1]
      else:
        aligned += "-" 

    seperated = aligned.split("-")
       
    score_gap = -1
    score_same = +1
    score_different = -1
    score = 0

    aligined2 = ""
    final = ""
    for i in aligned:
      if i == "":
        score += score_different
      elif i != "-" :
        final += "|"
        score += score_same
      
      else:
        final += "-"
        score += score_different

    print(seq11)
    print(final)
    print(seq2)

    return score

def using_bio_local(s1,s2):
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match = 1
    aligner.mismatch = -1
    aligner.gap_score = -1
    alignments = aligner.align(s1,s2)
    for alignment in alignments:
        print(alignment)
        print(alignment.score)

def using_bio_global(s1,s2):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match = 1
    aligner.mismatch = -1
    aligner.gap_score = -1
    alignments = aligner.align(s1,s2)
    for alignment in alignments:
        print(alignment)
        print(alignment.score)



#print(smith_waterman("GCCAT", "GCACT")) 


###############################################################################################################################################

      
def button_clicked():
    if var.get() == 1:
        
        res=needleman_wunsch(enter1.get().upper(), enter2.get().upper())
        global_alignment_brute(enter1.get().upper(), enter2.get().upper())
        print(greedy_align(enter1.get().upper(), enter2.get().upper()))
        using_bio_global(enter1.get().upper(), enter2.get().upper())
       
        label.config(text=res)
        list.append(res)
        listbox.insert(END, res)    

    elif var.get() == 2:
       res=smith_waterman(enter1.get().upper(), enter2.get().upper())
       using_bio_local(enter1.get().upper(), enter2.get().upper())
       label.config(text=res)
       list.append(res)
       listbox.insert(END, res)

       
###############################################################################################################################################

master = Tk()
master.title("Sequence aligment")
master.geometry("500x300")



bg = PhotoImage(file = "cell-and-molecular-biology.png")
# Show image using label
label1 = Label( master, image = bg)
label1.place(x = 0, y = 0)


#first label
Label(master,text="pairwise sequence aligment", font=("Arial", 12), fg="black", bg="#F0F0F8",borderwidth=2, relief="groove", padx=10, pady=10).place(x=145,y=10)

#first label
Label(master,text="Enter your sequences", font=("Arial", 12), fg="white", bg="#3232CD",borderwidth=2, relief="groove", padx=10, pady=10).place(x=165,y=60)

# Entry box for first sequence
Label(master,text="Enter your first sequence").place(x=182,y=120)

enter1= Entry(master)
enter1.place(x=188,y=145)

# Entry box for second sequence
Label(master,text="Enter your second sequence").place(x=174,y=170)
enter2= Entry(master)
enter2.place(x=185,y=195)


var = IntVar()

button1 = Radiobutton(master, text="global", variable=var, value=1, command=button_clicked)
button1.place(x=255,y=220)

button2 = Radiobutton(master, text="local", variable=var, value=2, command=button_clicked)
button2.place(x=180,y=220)

label = Label(master, text="")
label.place(x=195,y=250)

listbox = Listbox(master)


master.mainloop()