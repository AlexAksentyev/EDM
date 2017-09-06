function er = createrand(n)

ey=1;
er = randn(n,1)*ey/3;

for i=1:16
    if(er(i)>ey)
        er(i) = ey;
    end
    if(er(i)<-ey)
        er(i) = -ey;
    end
end

end