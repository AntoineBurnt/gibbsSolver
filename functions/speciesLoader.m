function [species,atoms] = speciesLoader(filename,reactants,source)

A = char(readlines(filename));
headers = [];

for i = 1:length(A(:,1))
    if contains(A(i,1),lettersPattern) || A(i,1) == '('
        headers = vertcat(headers,i);
    end
end

n_s = length(headers);

cards = cell(n_s,1);

for i = 1:n_s
    header = headers(i);
    if i == n_s
        next = length(A(:,1));
    else
        next = headers(i+1)-1;
    end
    cards{i} = A(header:next,:);
end

species = cellfun(@parse,cards);

toRemove = [species.phase] > 0;
for i = n_s:-1:1
    if toRemove(i)
        species(i) = [];
        n_s = n_s - 1;
    end
end


%Missing Textbook species
species(n_s+1) = struct('name',"C12H24",'intervals',[200 1000;1000 6000],'atoms',["C" "H"],'etta',[12 24],'phase',0,'weight',168.3226,'coefficients',[0 0 0.004991557	0.000100649	-1.59114E-08	-2.82088E-11	1.19647E-14	-19.543438	0.016904105;0 0 0.073484623	0	0	0	0	-47.675298	-0.36994817],'h_f',0,'h_f0',0);
species(n_s+2) = struct('name',"C6H12O",'intervals',[200 1000;1000 6000],'atoms',["C" "H" "O"],'etta',[6 12 1],'phase',0,'weight',100.1607,'coefficients',[0 0 0.46835	0.068027	-0.000039907	9.8998E-10	0	-36987	27.96;0 0 29.578	0	0	0	0	-45138.3	-124.789],'h_f',0,'h_f0',0);
n_s = n_s + 2;

if ~source
    species = load('data/textbook_species.mat').speciesCopy;
    n_s = length(species);
    
end   

species_dict = dictionary([species.name],1:n_s);
atoms = unique([species(species_dict(reactants)).atoms]);


for i = length(species):-1:1
    if any(~matches(species(i).atoms,atoms)) || species(i).phase > 0
        species(i) = [];
    end
end

end

function s = parse(card)

s.name = convertCharsToStrings(strip(card(1,1:16)));

%Parse temp intervals
n_intervals = (length(card(:,1))-2)/3;
s.intervals = zeros(n_intervals,2);

for i = 1:n_intervals
    s.intervals(i,:) = str2double(convertCharsToStrings(strsplit(strip(card(i*3,1:22)))));
end



%Parse chemical formula
formula_txt = strip(card(2,11:50));
s.atoms = transpose(convertCharsToStrings(extract(formula_txt,lettersPattern)));
ratios = str2double(convertCharsToStrings(strsplit(strip(erase(formula_txt,s.atoms)))));
s.etta = ratios(1:length(s.atoms));

s.phase = str2double(convertCharsToStrings(card(2,51:52)));

s.weight = str2double(convertCharsToStrings(card(2,53:65)));

s.h_f = str2double(convertCharsToStrings(card(2,66:80)));

s.h_f0 = str2double(convertCharsToStrings(card(3,66:80)));

coefficients = zeros(n_intervals,9);

for i = 1:n_intervals
    lines = card(1+ i*3:2 + i*3,:);

    for ii = 1:5
        coefficients(i,ii) = str2double(replace(convertCharsToStrings(lines(1,(ii-1)*16+1:ii*16)),'D','E'));
    end
    coefficients(i,6) = str2double(replace(convertCharsToStrings(lines(2,1:16)),'D','E'));
    coefficients(i,7) = str2double(replace(convertCharsToStrings(lines(2,17:32)),'D','E'));

    coefficients(i,8) = str2double(replace(convertCharsToStrings(lines(2,49:64)),'D','E'));
    coefficients(i,9) = str2double(replace(convertCharsToStrings(lines(2,65:80)),'D','E'));
end
s.coefficients = coefficients;
end