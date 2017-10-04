function solve_beam();

m = 4;
w1_0 = ones(m,1);
w2_0 = ones(m,1);
T_0 = ones(m,1);
pressure = ones(m+1,1);

dx = 1;
dt = 1;

G = @beam_model;
[y,term] = Broyden(G,w2_0)

    function g = beam_model(w2);

    % w1 - theta, (n-1)*1
    % w2 - d theta/d t, (n-1)*1
    % T - longitual force, (n-1)*1
    w1 = w1_0 + dt*w2;

    n = length(w1) + 1; % number of element

    % transformation matrix
    I = eye(n);
    I_addtop = I(:,2:end);
    I_adddown = I(:,1:end-1);
    % derivative matrix
    D = (- diag(ones(n,1),0) + diag(ones(n-1,1),1));
    D = D(1:end-1,:)./dx; % D is (n-1)*n

    dw1ds = (I_adddown * D) * (I_addtop * w1);
    dw1ds_2 = (I_adddown * D) * dw1ds;
    dw1ds_3 = zeros(n,1);
    dw1ds_3(2:end) = D * dw1ds_2;

    F = @inextensible_fn;
    Tn = Broyden(F, T_0);
    Tn = I_adddown * Tn;

    dw1ds_3(1) = Tn(1) * dw1ds(1) - pressure(1);
    dw1ds_4 = D * dw1ds_3; % (n-1)*1 
    dpds = D * pressure; % (n-1)*1 
    dTds = D * Tn; % (n-1)*1 

    dw1dt = w2;
    dw2dt = -dpds - dw1ds_4 + (Tn(2:end) ...
        + dw1ds(2:end).*dw1ds(2:end)).*dw1ds_2(2:end)...
        + 2 * dTds .* dw1ds(2:end);
    g = -w2 + w2_0 + dt * dw2dt;

        function f = inextensible_fn(Tx);
            dw1ds_3(1) = Tx(1) * dw1ds(1) - pressure(1);

            dTds = zeros(n,1);
            dTds(2:end) = D * I_adddown * Tx;
            dTds(1) = - dw1ds(1) * dw1ds_2(1);
            dTds_2 = D * dTds;

            f = dTds_2 - Tx.*dw1ds(1:end-1).*dw1ds(1:end-1) ...
                + pressure(1:end-1).*dw1ds(1:end-1) ...
                + 2*dw1ds(1:end-1).*dw1ds_3(1:end-1) ...
                + dw1ds_2(1:end-1).*dw1ds_2(1:end-1) ...
                + w2.*w2;
        end
    end
end
