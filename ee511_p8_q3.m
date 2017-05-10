%%--------------------------------------------------------------------------
%%Project-8:: Question - 3
%%To simulate simulated annealing
%%Author                Date               Revision
%%Rajasekar Raja     05/07/17         Initial Revision
%%--------------------------------------------------------------------------
function [] = ee511_p8_q3(no_of_samples)
min_limit = -500;
max_limit = 500;
equ_spaced_pts1 = linspace(min_limit,max_limit); %Row vector of 100 linearly equally spaced points between -512 and 512
equ_spaced_pts2 = equ_spaced_pts1;
[equ_spaced_pts1,equ_spaced_pts2] = meshgrid(equ_spaced_pts1,equ_spaced_pts2);
func = 418.9829*2 - equ_spaced_pts1.*sin(sqrt(abs(equ_spaced_pts1))) - equ_spaced_pts2.*sin(sqrt(abs(equ_spaced_pts2)));
figure(1);
contour(equ_spaced_pts1,equ_spaced_pts2,func);
colorbar;

figure(2)
mesh(equ_spaced_pts1,equ_spaced_pts2,func);
colorbar;
xlim([-500 500]);
ylim([-500 500]);

matrix_x=zeros(no_of_samples);
matrix_y=zeros(no_of_samples);
matrix_Z=zeros(no_of_samples);
value1=zeros(no_of_samples);
value2=zeros(no_of_samples);
matrix_Z=100;
for iter = 1 :no_of_samples
    value1(iter+1) = matrix_x(iter) +normrnd(0,10);
    value2(iter+1) =matrix_y(iter) + normrnd(0,10);
    val1 = 418.9829*2 - matrix_x(iter)*sin(sqrt(abs(matrix_x(iter)))) - matrix_y(iter)*sin(sqrt(abs(matrix_y(iter))));
    val2 = 418.9829*2 - value1(iter+1)*sin(sqrt(abs(value1(iter+1)))) - value2(iter+1)*sin(sqrt(abs(value2(iter+1))));
    alpha = exp(val1 -val2)/ matrix_Z(iter);
    if ((val2 <= val1) || (rand(1)<alpha))
        matrix_x(iter+1) = value1(iter+1);
        matrix_y(iter+1) = value2(iter+1);     
    else
        matrix_x(iter+1)=  matrix_x(iter);
        matrix_y(iter+1) = matrix_y(iter); 
    end
    matrix_Z(iter+1) = 100/log(iter+1);
end
value1(no_of_samples);
value2(no_of_samples);
end