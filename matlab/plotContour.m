A=importdata('..\VScode\ReqSeqence.txt');
[rows, cols] = size(A);
totalRows = rows;
%%
%ȥ���쳣ֵ
ave = mean(A);
variance = var(A);
posLimit = ave + 3*variance;
negLimit = ave - 3*variance;
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

%S = sum(week,2);
[rows, cols] = size(week);
timeWeek = linspace(1, rows, rows);

figure
for i = 1:15
    subplot(5,3,i);
    plot(timeWeek',week(:,i),'*-');
end

figure
for i = 1:15
    subplot(5,3,i);
    plot(time',A(:,i),'*-');
end

figure
for i = 1:15
    subplot(5,3,i);
    plot(time',B(:,i),'*-');
end

figure 
for i=2:totalRows
    S(i) = S(i) + S(i -1);
end
plot(time', S, '.-');


































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
   
