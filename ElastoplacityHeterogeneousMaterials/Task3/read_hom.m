A=dlmread('task2.out.hom_rand');
t=A(:,1);
epsilon=A(:,2:7);
sigma=A(:,8:13);

rsigma_v  =von_Mises(sigma  (:,1),sigma  (:,2),sigma  (:,3),sigma  (:,4),sigma  (:,5),sigma  (:,6));
repsilon_v=von_Mises(epsilon(:,1),epsilon(:,2),epsilon(:,3),epsilon(:,4),epsilon(:,5),epsilon(:,6));