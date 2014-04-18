function [energy, inds, ngbr] = get_energy_nab(alpha, beta, labels, unary_pot, vertC, horC, metric)
% This function calculates energy for given label configuration
% INPUT:
%    labels: N-by-M array of integer, labels for pixels
%    unary_pot: N-by-M-by-K array of double, unary potentials
%    pair_pot: N-by-M-by-K array of double, pixel brightness
%
% OUTPUT:
%    energy: double

    inds = [];
    ngbr = [];
    [N, M] = size(labels);
    energy = 0;
    for i = 1 : N
        for j = 1 : M
            if labels(i, j) ~= alpha && labels(i, j) ~= beta
                energy = energy + unary_pot(i, j, labels(i, j));
                inds = [inds, i + N * (j - 1)]; %#ok
                if (i < N)
                    if labels(i + 1, j) ~= alpha && labels(i + 1, j) ~= beta
                        energy = energy + vertC(i, j) * metric(labels(i, j), ...
                            labels(i + 1, j));
                    else
                        ngbr = [ngbr, i + N * (j - 1)]; %#ok
                    end
                end
                if (j < M)
                    if labels(i, j + 1) ~= alpha && labels(i, j + 1) ~= beta
                        energy = energy + horC(i, j) * metric(labels(i, j), ...
                            labels(i, j + 1));
                    else
                        ngbr = [ngbr, i + N * (j - 1)]; %#ok
                    end
                end
            else
                if (i < N)
                    if labels(i + 1, j) ~= alpha && labels(i + 1, j) ~= beta
                        ngbr = [ngbr, i + 1 + N * (j - 1)]; %#ok
                    end
                end
                if (j < M)
                    if labels(i, j + 1) ~= alpha && labels(i, j + 1) ~= beta
                        ngbr = [ngbr, i + N * j]; %#ok
                    end
                end
            end
        end
    end
end