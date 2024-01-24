function S1Img = S1addOutput(cI, S1res, nFaces, seeMadeImg)

clear matPerFreq1 matPerFreq

nStimPerDir = length(cI);

scrsz = get(0,'ScreenSize');
fsz = 1.1;
duration_pause = 0.7; % can be shortened or prologed at will


for j = 1:nStimPerDir % nStim per directory

    cnt = 0;
    for k = 1:8 % nBands

        for l = 1:4 % oriented filters

            cnt = cnt+1;

            SeeStim1 = S1res{j,k}{1}{1,l};
            SeeStim2 = S1res{j,k}{2}{1,l};

            if seeMadeImg
                display(['j = ' num2str(j)])

                figure('Position',[1 scrsz(4)/fsz scrsz(3)/fsz scrsz(4)/fsz])

                subplot(1,2,1),imagesc(SeeStim1),axis equal tight
                subplot(1,2,2),imagesc(SeeStim2),axis equal tight

                pause(duration_pause), close
            end


            matPerFreq1(:,:,cnt) = SeeStim1;
            matPerFreq2(:,:,cnt) = SeeStim2;

        end

    end

    SumFreq{j} = sum(matPerFreq1,3) + sum(matPerFreq2,3);
end
S1Img = reshape(SumFreq,[5, nFaces]);
S1Img = S1Img';
