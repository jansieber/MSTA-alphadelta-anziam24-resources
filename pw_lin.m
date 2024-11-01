function mat=pw_lin(tc,b)
mat=exp(-b*tc)*(eye(2)+[b*tc,tc; -b^2*tc,-b*tc]);
end
