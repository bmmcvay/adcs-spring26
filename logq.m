function mat=logq(q)
    theta=acos(q(1));
    r=q(2:4)/norm(q(2:4));
    mat=theta*r;
end