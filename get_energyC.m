function energy = get_energyC(labels, unary_pot, vertC, horC, metric)
% This function calculates energy for given label configuration
% INPUT:
%    labels: N-by-M array of integer, labels for pixels
%    unary_pot: N-by-M-by-K array of double, unary potentials
%    pair_pot: N-by-M-by-K array of double, pixel brightness
%
% OUTPUT:
%    energy: double

    [N, M] = size(labels);
    energy = 0;
    for i = 1 : N
        for j = 1 : M
            energy = energy + unary_pot(i, j, labels(i, j));
            if (i < N)
                energy = energy + vertC(i, j) * metric(labels(i, j), ...
                    labels(i + 1, j));
            end
            if (j < M)
                energy = energy + horC(i, j) * metric(labels(i, j), ...
                    labels(i, j + 1));
            end
        end
    end
end