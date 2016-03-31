
% U = {rand(5,3), rand(4,2), rand(3,1)};
% core1 = sptenrand([3 2 1],3);
% Y = ttensor(core1,U);
% Y = full(Y);
% icdlfullpath = mfilename('fullpath');
% function [] = test_full(d)
% C = mfilename('class') ;
% end
fid = fopen('test', 'wb');
a = [1,0,0,1,1,0,1,1,0,1,1,1,0,0,1,1];
a = uint8(a);
fwrite(fid, a(1:8), 'ubit8');
fwrite(fid, a(9:16), 'ubit8');
fclose(fid);