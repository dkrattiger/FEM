function [K,M] = master_mass_stiffness_3DE_BlochOperator(kappa,Mesh_info)

% Recover Mesh info
Ds              = Mesh_info.Ds;
rhos            = Mesh_info.rhos;
pattern_vec     = Mesh_info.pattern_vec;
n               = Mesh_info.n;
n_dof           = Mesh_info.n_dof;
elenodes        = Mesh_info.elenodes;
coordinates     = Mesh_info.coordinates;

n_kap  = size(kappa,2);
n_eles = size(elenodes,1);

i_col = zeros((3*n^3)^2*n_eles,1);
i_row = zeros((3*n^3)^2*n_eles,1);
K_val = zeros((3*n^3)^2*n_eles,1);
M_val = zeros((3*n^3)^2*n_eles,1);

for i = 1:n_eles
    i
    dofs_in_ele = zeros(3*n^3,1);
    for j = 1:length(elenodes(i,:))
        en = elenodes(i,j);
        dofs_in_ele(3*(j-1)+1:3*j) = [3*(en-1)+1:3*en];
    end

    rho = rhos{pattern_vec(i)};
    D   = Ds{pattern_vec(i)};

    [K_ele,M_ele] = element_mass_stiffness_3DE(D,rho,n,kappa,...
        coordinates(elenodes(i,:),:));
    
    
    i_row_ele = repmat(dofs_in_ele,[1,3*n^3]);
    i_col_ele = i_row_ele';
    
    i_row_ele = i_row_ele(:);
    i_col_ele = i_col_ele(:);
    
    i_col(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = i_col_ele;
    i_row(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = i_row_ele;
    
    K_val(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = K_ele(:);
    M_val(((3*n^3)^2*(i-1)+1):((3*n^3)^2*i)) = M_ele(:);

end

K = sparse(i_row,i_col,K_val);
M = sparse(i_row,i_col,M_val);

K = (1/2)*(K+K');
M = (1/2)*(M+M');