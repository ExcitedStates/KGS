%Import and convert Kinari hbond files

clear all; close all; clc;

proteinNameA='3K0N_A_H92F';
proteinNameB='3K0N_B_H92F';
comboName='3K0N_AB_H92';
multiChain = 0;

%A=importdata(['bondFiles/',proteinNameA,'.hBonds.bnd.knr']);
A=importdata(['bondFiles/',proteinNameA,'.HBonds.Pruned.bnd.knr']);
% A=importdata(['bondFiles/',proteinName1]);

% proteinNameA='gas-b2ar-inac-new';
% proteinNameB='gas-b2ar-acti-new';
% comboName='gas-b2ar';
% A=importdata(['bondFiles/hbonds_gas_b2ar_inac.in']);


% KINARI files
atom1=A(:,2);
atom2=A(:,4);
% energy=A(:,5);

%Other files
% atom1=A(:,1);
% atom2=A(:,2);
% energy=A(:,3);

matrixA=[atom1  atom2  ] ; %energy];


%load file
B=importdata(['bondFiles/',proteinNameB,'.HBonds.Pruned.bnd.knr']);
% B=importdata(['bondFiles/3K0N_AB_H92F_hb2_hBond_pairs.txt']);
% B=importdata(['bondFiles/',proteinName2]);
% B=importdata(['bondFiles/hbonds_gas_b2ar_acti.in']);

% KINARI files
atom1=B(:,2);
atom2=B(:,4);
% energy=B(:,5);

%Other files
% atom1=B(:,1);
% atom2=B(:,2);
% energy=B(:,3);

matrixB=[atom1  atom2 ]; %  energy];

%% Import PDB files

pdbA = pdbread(['bondFiles/',proteinNameA,'.Processed.pdb.knr']);
%pdbA = pdbread(['bondFiles/',proteinNameA,'_Processed.pdb']);
% pdbA = pdbread(['bondFiles/',proteinNameA,'.pdb']);
pdbA = pdbA.Model.Atom;

pdbB = pdbread(['bondFiles/',proteinNameB,'.Processed.pdb.knr']);
%pdbB = pdbread(['bondFiles/',proteinNameB,'_Processed.pdb']);
% pdbB = pdbread(['bondFiles/',proteinNameB,'.pdb']);
pdbB = pdbB.Model.Atom;

minLength=min(length(matrixA),length(matrixB));
combo=zeros(minLength,2); %3);

%% Sort hbonds (first id < second id)
for i=1:length(matrixA)
    if(matrixA(i,2) < matrixA(i,1) )
        tmp=matrixA(i,2);
        matrixA(i,2)=matrixA(i,1);
        matrixA(i,1)=tmp;
    end
end

for i=1:length(matrixB)
    if(matrixB(i,2) < matrixB(i,1) )
        tmp=matrixB(i,2);
        matrixB(i,2)=matrixB(i,1);
        matrixB(i,1)=tmp;
    end
end

count=0;

for i=1:length(matrixA)
    id1A = matrixA(i,1);
    id2A = matrixA(i,2);
    ind1A = find([pdbA.AtomSerNo] == id1A);
    ind2A = find([pdbA.AtomSerNo] == id2A);

    if ~isempty(ind1A) && ~isempty(ind2A)
        atom1A = pdbA(1,ind1A);
        atom2A = pdbA(1,ind2A);

        resId1A = atom1A.resSeq;
        atomName1A = atom1A.AtomName;
        chainId1 = atom1A.chainID;
        resId2A = atom2A.resSeq;
        atomName2A = atom2A.AtomName;
        chainId2 = atom2A.chainID;

        resSel1B = ([pdbB.resSeq] == resId1A);
        resSel1B = pdbB(1,resSel1B);
        resSel2B = ([pdbB.resSeq] == resId2A);
        resSel2B = pdbB(1,resSel2B);

        found1=-1;
        found2=-1;

        for j=1:length(resSel1B)
            if(resSel1B(1,j).chainID == chainId1 && strcmp(resSel1B(1,j).AtomName,atomName1A))
                id1B = resSel1B(1,j).AtomSerNo;
                found1=1;
                break;
            end
        end
        for j=1:length(resSel2B)
            if(resSel2B(1,j).chainID == chainId2 && strcmp(resSel2B(1,j).AtomName,atomName2A))
                id2B = resSel2B(1,j).AtomSerNo;
                found2=1;
                break;
            end
        end

        if(found1 == 1 && found2 == 1)
            atom1B = find([pdbB.AtomSerNo] == id1B);     
            atom2B = find([pdbB.AtomSerNo] == id2B);

            ind1B = find(matrixB(:,1)==id1B);
            for finalInd = 1:length(ind1B)
                foundB = matrixB(ind1B(finalInd),2)== id2B;  
                if(foundB)
                    count=count+1;
%                     combo(count,:)=[id1A, id2A, matrixA(i,3)];
                    combo(count,:)=[id1A, id2A];
                end
            end
        end
    end  
end

output=combo(1:count,:)';
disp(['Found ',num2str(count),' common hbonds.']);
fid=fopen(['bondFiles/',comboName,'_hBond_pairs.txt'],'w');

% fprintf(fid, '%4d %4d %12.8f\n',output);
fprintf(fid, '%4d %4d \n',output);

fclose(fid);


%% Identify inter-chain hydrogen bonds

if multiChain == 1
    
output=output';

chainBorderA = 4818; %5755; %
chainBorderB = 4818; %5755; %

% chainBorder2 = 6677;
% 
% interval = find(output(:,1)>=chainBorder1 & output(:,1)<=chainBorder2);
% check = find(output(interval,2) <=chainBorder1 | output(interval,2) >=chainBorder2);
% disp(output(interval(check),:))

interval = find(output(:,1)<=chainBorderA);
check = find(output(interval,2) >chainBorderA);
disp('Common inter-chain bonds:');
disp(output(interval(check),1:2))

%% In Matrix A 
% interval = find(matrixA(:,1)<=chainBorderA);
% check = find(matrixA(interval,2) >chainBorderA);
% disp('inter-chain bonds in matrix A:');
% disp(matrixA(interval(check),1:2))
% 
% %% In Matrix B
% interval = find(matrixB(:,1)<=chainBorderB);
% check = find(matrixB(interval,2) >chainBorderB);
% disp('inter-chain bonds in matrix B:');
% disp(matrixB(interval(check),1:2))

end
