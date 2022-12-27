%%
k = 290;
num = 0;
for  i1= 0:k
    for i2 =0:k
        for i3 =0:k
            for i4 = 0:k
                if i1+2*i2+3*i3+4*i4 == k
                    num = num +1;
                end
            end
        end
    end
end