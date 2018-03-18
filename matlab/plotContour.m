A=importdata('..\VScode\ReqSeqence.txt');
%A=importdata('..\VScode\data\TestData_2015.2.20_2015.2.27.txt');
[rows, cols] = size(A);
totalRows = rows;
%%
%È¥³ýÒì³£Öµ
ave = mean(A);
variance = var(A);
posLimit = ave + 2*variance;
negLimit = ave - 2*variance;
for r = 1:rows
    for c = 1:cols
        if A(r,c) > posLimit(c) || A(r,c) < negLimit(c) 
            if r~=1 && r~= rows
                A(r,c) = (A(r - 1,c) + A(r + 1, c))/2;
            end
        end
    end
end


S = sum(A,2);
%%
B = A;
for i = 2 : rows
    B(i,:) = B(i,:) + B(i-1,:);
end
time = linspace(1, rows, rows);
week = [];
interval = 7;
for i = 1:interval:rows - interval
    sums = 0;
   for j = 1:interval
       sums = sums + A(j,:);       
   end     
   week = [week; sum(A(i:i + interval - 1,:))];
end
%%

%S = sum(week,2);
[rows, cols] = size(week);
for i=2:rows
    week(i,:) = week(i,:) + week(i - 1,:);
end

timeWeek = linspace(1, rows, rows);

figure(1)
for i = 1:15
    subplot(5,3,i);
    plot(timeWeek',week(:,i),'*-');
end

% figure
% for i = 1:15
%     subplot(5,3,i);
%     plot(time',A(:,i),'*-');
% end
% 
% figure
% for i = 1:15
%     subplot(5,3,i);
%     plot(time',B(:,i),'*-');
% end

% figure 
% for i=2:totalRows
%     S(i) = S(i) + S(i -1);
% end
% plot(time', S, '.-');
coef = [];
seqX = linspace(1, totalRows, totalRows);
for i = 1 : cols
    p = polyfit(seqX', A(:,i), 3);
    coef = [coef;p];
end
for day = 43:49
    for i = 1 : cols
        predict = polyval(coef(i,:), day);
        A(day,i) = predict;
    end
end

figure(2)
timeSeq = linspace(1, totalRows + 7, totalRows + 7);
for i = 1:15
    subplot(5,3,i);
    plot(timeSeq',A(:,i),'*-');
end
predVal = sum(A(43:end,:),1)

































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
   
