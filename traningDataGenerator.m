%%Continuously run the FinalProject script to generate traning data

tic;
traningSetLength = 1000;

for i=1:traningSetLength
    fprintf("Generating traning data --- <strong>%.1f%%</strong> complete\n", 100*(i-1)/(traningSetLength-1));
    FinalProject
end
toc;