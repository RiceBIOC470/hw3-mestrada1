%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

% See repository for snapshot of solution

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

accession1 = 'NM_002746';
erk1_data = getgenbank(accession1);
erk1seq = erk1_data.Sequence; %1902 base pairs

accession2 = 'NM_002745';
erk2_data = getgenbank(accession2);
erk2seq = erk2_data.Sequence;

[score, align, start] = swalign(erk1seq, erk2seq, 'Alphabet', 'nt', 'Showscore', true);

base_matches = count(align(2,:), '|');
frac_base_pairs = base_matches/(length(erk1seq));

% Fraction of base pairs in ERK1 that align to ERK2 is 0.5536

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

accession1pro = erk1_data.CDS.protein_id;
accession2pro = erk2_data.CDS.protein_id;
erk1pro_data = getgenpept(accession1pro);
erk2pro_data = getgenpept(accession2pro);
erk1proseq = erk1pro_data.Sequence; %379 aa
erk2proseq = erk2pro_data.Sequence; %360 aa

[score, align, start] = swalign(erk1proseq, erk2proseq, 'Alphabet', 'AA', 'Showscore', true);

aa_matches = count(align(2,:), '|');
frac_aa = aa_matches/(length(erk1proseq));

% Fraction of amino acids in ERK1 that align to ERK2 is 0.8057

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

% DNA SEQUENCE ALIGNMENT

% Accession # for mus musculus ERK1: NM_011952.2 
% Accession # for mus musculus ERK2: NM_011949.3 

accession1mus = 'NM_11952.2';
accession2mus = 'NM_011949.3';
erk1mus_data = getgenbank(accession1mus);
erk2mus_data = getgenbank(accession2mus);
erk1musseq = erk1mus_data.Sequence;
erk2musseq = erk2mus_data.Sequence;

% ERK1
[score, align, start] = swalign(erk1seq, erk1musseq, 'Alphabet', 'nt', 'Showscore', true);
base_matches_erk1 = count(align(2,:), '|');
frac1_erk1 = base_matches_erk1/(length(erk1seq));
% Fraction of ERK1 human base pairs that align with ERK1 mouse base pairs:
% 0.7939
frac2_erk1 = base_matches_erk1/(length(erk1musseq));
% Fraction of ERK1 mouse base pairs that align with ERK1 human base pairs:
% 0.8521
% These sequences are pretty similar

% ERK2
[score, align, start] = swalign(erk2seq, erk2musseq, 'Alphabet', 'nt', 'Showscore', true);
base_matches_erk2 = count(align(2,:), '|');
frac1_erk2 = base_matches_erk2/(length(erk2seq));
% Fraction of ERK2 human base pairs that align with ERK2 mouse base pairs:
% 0.6611
frac2_erk2 = base_matches_erk2/(length(erk2musseq));
% Fraction of ERK2 mouse base pairs that align with ERK2 human base pairs:
% 0.7670
% These sequences are decently similar (not as similar as the ERK1
% sequences)

% PROTEIN SEQUENCE ALIGNMENT

accession1muspro = erk1mus_data.CDS.protein_id;
accession2muspro = erk2mus_data.CDS.protein_id;
erk1mus_prodata = getgenpept(accession1muspro);
erk2mus_prodata = getgenpept(accession2muspro);
erk1mus_proseq = erk1mus_prodata.Sequence;
erk2mus_proseq = erk2mus_prodata.Sequence;

% ERK1 
[score, align, start] = swalign(erk1mus_proseq, erk1proseq, 'Alphabet', 'AA', 'Showscore', true);
aa_match_erk1 = count(align(2,:),'|');
frac1_aa_erk1 = aa_match_erk1/(length(erk1proseq));
% Fraction of ERK1 human amino acids that align with ERK1 mouse amino
% acids: 0.9683
frac2_aa_erk1 = aa_match_erk1/(length(erk1mus_proseq));
% Fraction of ERK1 mouse amino acids that align with ERK1 human amino
% acids: 0.9658
% These amino acid sequences are very similar

% ERK2 
[score, align, start] = swalign(erk2proseq, erk2mus_proseq, 'Alphabet', 'AA', 'Showscore', true);
aa_match_erk2 = count(align(2,:), '|');
frac1_aa_erk2 = aa_match_erk2/(length(erk2proseq));
% Fraction of ERK2 human amino acids that align with ERK2 mouse amino
% acids: 0.9861
frac2_aa_erk2 = aa_match_erk2/(length(erk2mus_proseq));
% Fraction of ERK2 mouse amino acids that align with ERK2 human amino
% acids: 0.9916
% These amino acid sequences are very similar 

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

% See blasthits.m file in repository

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

% See blast_hscore_matches.m file in repository

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

% Human genome gene: Hypoxia inducible factor 1 alpha subunit (HIF1A);
% Accession #: NM_001530 

human_accession = 'NM_001530';
[human_match, nonhuman_match] = blast_hscore_matches(human_accession);
human_match
human_match =

    'NM_001243084'

nonhuman_match
nonhuman_match =

    'XM_003831627'
    
humanmatch_data = getgenbank(human_match);
humanmatch_data.Definition
ans =

    'Homo sapiens hypoxia inducible factor 1 alpha subunit (HIF1A), transcript variant 3, mRNA.'
% The human match makes sense; HIF1A aligns with HIF1A, mRNA transcript variant 3    

nonhumanmatch_data = getgenbank(nonhuman_match);
nonhumanmatch_data.Definition
ans =

    'PREDICTED: Pan paniscus hypoxia inducible factor 1, alpha subunit (basic helix-loop-helix transcription factor) (HIF1A), transcript variant X2, mRNA.'
% The nonhuman match makes sense; this gene comes from Pan paniscus and is
% also a hypoxia inducible factor 1, alpha subunit gene

% Nonhuman genome (Mus musculus) gene: Kirsten rat sarcoma viral oncogene
% homolog (KRAS); Accession #: NM_021284

kras_accession = 'NM_021284';
[human_match, nonhuman_match] = blast_hscore_matches(kras_accession);
human_match
human_match

    'No close high scoring match in human DNA/RNA'
% The lack of a high-scoring human match makes sense; the KRAS gene has 
% homologs in many different species, all of which have higher scoring
% sequence alignments with Mus musculus KRAS than any human gene has with
% KRAS.

nonhuman_match
nonhuman_match =

    'XM_006506918'
nonhumanmatch_data2 = getgenbank(nonhuman_match);
nonhumanmatch_data2.Definition

ans =

    'PREDICTED: Mus musculus Kirsten rat sarcoma viral oncogene homolog (Kras), transcript variant X1, mRNA.'
% The nonhuman match makes sense; this match is another Mus musculus KRAS
% gene, but it is slightly modified as it is a transcript variant. 