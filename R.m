function mat=R(q)
    mat=[q(1) -q(2:4)';q(2:4) q(1)*eye(3)-hat(q(2:4))];
end