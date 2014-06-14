function Estimate(ChrNameList, GTFPath, OutputFilePath)

[ChrNames] = textread(ChrNameList, '%s');
ChrNum = length(ChrNames);

for chrLoopCnt = 1:ChrNum
    ChrNames{chrLoopCnt}
    filepath = sprintf('%s%s/GeneList.txt', GTFPath, ChrNames{chrLoopCnt});
    [GeneNames, ensemblIden, symbol, transNum] = textread(filepath, '%s%s%s%n');
    GeneNum = length(GeneNames);
    for geneLoopCnt = 1:GeneNum
        outputfilename = sprintf('%s%s/temp/result_%s.txt', OutputFilePath, ChrNames{chrLoopCnt}, GeneNames{geneLoopCnt});
        file = fopen(outputfilename,'w+t');
        TransFrag = sprintf('%s%s/temp/structure_%s.txt', OutputFilePath, ChrNames{chrLoopCnt}, GeneNames{geneLoopCnt});
        [M]=textread(TransFrag);
        [m,n] = size(M);
        M = M(:, 1:(n-1));
        Coverage = sprintf('%s%s/temp/coverage_%s.txt', OutputFilePath, ChrNames{chrLoopCnt}, GeneNames{geneLoopCnt});
        c = textread(Coverage);        
        T = lsqnonneg(M,c);        
        for i = 1:length(T)
            fprintf(file, '%f\n', T(i));
        end;
        fclose(file);
    end;
end;








