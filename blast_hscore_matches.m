function [human_match, nonhuman_match] = blast_hscore_matches(accession)
blast_hits_accessions = blasthits(accession);
human_match = [];
for ii=2:50 %Exclude first hit (usually alignment with own sequence)
    current_access = blast_hits_accessions{ii};
    data = getgenbank(current_access);
    if strfind(data.Source, 'Homo sapiens') == 1
        human_match = data.Accession;
        break
    end
end
if isempty(human_match) == 1 
    human_match = 'No close high scoring match in human DNA/RNA';
end
for ii=2:50
    current_access = blast_hits_accessions{ii};
    data = getgenbank(current_access);
    if isempty(strfind(data.Source, 'Homo sapiens')) == 1
        nonhuman_match = data.Accession;
        break
    end
end
end