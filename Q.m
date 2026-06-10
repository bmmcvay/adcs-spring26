function mat=Q(q)
    mat=H'*L(q)*R(q)'*H();
end