# Prediction-of-polypeptide-hydrolysis
The use of mass spectrometry to verify the structure of unknown peptides is often done using a 'shotgun approach' whereby the peptides are hydrolyzed randomly into fragments by proteases and then measure the sample by using a high-resolution mass spectrometry. However, complex peptides especially those with ring structure, have a large number of hydrolysates, it is a labor-intensive task to determine whether the molecular ions given by mass spectrometry belong to the target peptides. This program provides the ability to predict all possible hydrolysates of a given peptide and order them according to precise molecular weight.

Replace ‘testnode’ and ‘testbond’ list in the main.py which include amino acids and bonds information. Input amino acids in the ‘testnode’ list in order and input the list of bonds consists of a ternary array. The first and second number of array represents the sequence number of the two amino acids to which the bond is attached and the third number represents the type of bond. Then you can run the program. In the current version, only peptide bonds (type 1) and disulfide bonds (type 2) are supported. The output of the program uses IUPAC condense format to represent the structure of the polypeptide.

The principle of this program is to treat the amino acid monomer as a node of a graph, and the bond between amino acids as an edge of the graph. In this way, the problem of predicting all hydrolyzed fragments of a polypeptide is equivalent to the problem of enumerating all the subgraphs of the above graph. In addition, the algorithm used in this program is optimized by the characteristics of peptides. Since there are usually many chain structures in a polypeptide, each single chain is treated as a supernode. First, we enumerate all the subgraphs of the graph composed of supernodes. Then, we find the supernodes connected to these subgraphs which will be treated as incomplete chains. The final result is obtained by sequentially chopping the bonds on these single chains and merging them with the subgraph. This method improves the performance of the program by avoiding the need to store all graph structures in memory.

There is currently no license for this program. You can use it for calculations, but cannot disseminate it in any way or use it for commercial services. If you have any need, please contact me by sending email to xuezixuan01@gmail.com 

The C++ version of this program optimizes the algorithm of exhaustive subgraphs and some data structures, which makes it faster and supports the generation of tens of millions of fragments. However, its ownership belongs to Che-Hung Micro Technology (Shanghai) Co., Ltd.(hzwtech)
 and it will be a feature of hzwtech 's subsequent products. If necessary, please pay attention to HZWtech 's product release In the future.
