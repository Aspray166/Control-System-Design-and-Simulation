function zk(a,b)
num=a;den=a+b;
g=tf(num,den);
[z,p,k]=tf2zp(num,den);
g2=zpk(z,p,k);
g3=tf2ss(num,den);
[num,den]=tfdata(g,'v');
[r,p,k]=residue(num,den);
end
a=[0,0,10]
b=conv([1,0],[1,1])