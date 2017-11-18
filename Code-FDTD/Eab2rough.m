function Q=Eab2rough(E_threshold,E_ab,Ez_ab,roughness_function) %#ok<INUSL>
E_ab=1e-6*E_ab;
Q=roughness_function;
Q(E_ab>=E_threshold)=0;
end