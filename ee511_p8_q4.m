%%--------------------------------------------------------------------------
%%Project-8:: Question - 4
%%To simulate Travelling salesman problem using Simulated Annealing
%%Author                Date               Revision
%%Rajasekar Raja     05/07/17         Initial Revision
%%--------------------------------------------------------------------------
clear all;
clc;
content = fileread('uscap_xy.txt'); 
data = textscan( content, '%f %f%*[^\n]','HeaderLines', 0) ;
x = data{1};
y = data{2};
n= 48;
city = [x y];

distance = pdist2(city, city);
num_iter = 100000;
c = 100;
p = [4:n 1:3];
len = 0;
for a1 = 1:n-1
    len = len + distance(p(a1),p(a1+1));
end
len = len + distance(p(n),p(1));
pathHistory = zeros(num_iter,n);
lenHistory = zeros(1,n);

% Plotting intial path
figure(1)
plot(city(:,1), city(:,2), 'ro');
xlim([min(x)-1 max(x)+1]);
ylim([min(y)-1 max(y)+1]);
hold on
line(city([p(:); p(1)],1), city([p(:); p(1)],2));
title('Initial path');
hold off

count = 0;
while(count<num_iter)
    count = count + 1;
    % Create path p2 by randomly swap two cities
    swap_index = randsample(n,2);
    p2 = p;
    temp = p2(swap_index(1));
    p2(swap_index(1)) = p2(swap_index(2));
    p2(swap_index(2)) = temp;
    % Cost of p2
    len2 = 0;
    for a1 = 1:n-1
        len2 = len2 + distance(p2(a1),p2(a1+1));
    end
    len2 = len2 + distance(p2(n),p2(1));
    q = (1+count)^((len - len2)/c);
    if len2 - len <= 0
        p = p2;
        len = len2;
    else
        if rand <= q
            p = p2;
            len = len2;
        end
    end
    pathHistory(count,:) = p;
    lenHistory(count) = len;
end

figure(2)
plot(1:num_iter, lenHistory, 'linewidth',2);
title('Length of path in each iteration')

figure(3)
plot(city(:,1), city(:,2), 'ro');
xlim([min(x)-1 max(x)+1]);
ylim([min(y)-1 max(y)+1]);
hold on
line(city([p(:); p(1)],1), city([p(:); p(1)],2));
title('Final path of Travelling Salesman Problem using Simulated Annealing');
hold off

fprintf('\nDistance of the shortest path that visits each of the 48 state capital cities is: %f', len);