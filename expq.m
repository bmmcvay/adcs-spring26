function mat=expq(phi)
    mat=[cos(norm(phi));phi*sinc(norm(phi)/pi)];
end