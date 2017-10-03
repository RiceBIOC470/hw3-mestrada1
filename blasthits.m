function output = blasthits(accession_num, N)

% NCBI identifies nucleotide sequences with a two-letter alphabetical
% prefix and protein sequences with a three-letter alphabetical prefix;
% here, I utilize this distinction to apply this function to both mRNA and
% protein accession numbers.

if exist('N') == 0
    N = 50; % defaulted to 50 hits if input not specified
end
if isstrprop(accession_num(3), 'alpha') == 0
    gene_data = getgenbank(accession_num);
    gene_sequence = gene_data.Sequence;
    [requestID, requestTime] = blastncbi(gene_sequence, 'blastn', 'Database', 'refseq_rna');
    blast_data = getblast(requestID, 'WaitTime', requestTime);
else 
    protein_data = getgenpept(accession_num);
    pro_sequence = protein_data.Sequence;
    [requestID, requestTime] = blastncbi(pro_sequence, 'blastp', 'Database', 'refseq_protein');
    blast_data = getblast(requestID, 'WaitTime', requestTime);
end
cell_array{1,N} = [];
for ii = 1:N
    blast_hit = blast_data.Hits(ii).Name;
    indeces = strfind(blast_hit, '|');
    accession = blast_data.Hits(ii).Name((indeces(3)+1):(indeces(4)-1));
    cell_array{ii} = accession;
end
output = cell_array;
