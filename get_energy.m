function energy = get_energy(labels, unary_pot, intensity) 
            %vertC, horC, metric)
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
                energy = energy + ...
                    get_pp(i + N * (j - 1), i + 1 + N * (j - 1), ...
                    labels(i, j), labels(i + 1, j), intensity);
                    abs(intensity(i, j, labels(i, j)) - intensity(i, j, labels(i + 1, j))) + ...
                    abs(intensity(i + 1, j, labels(i, j)) - intensity(i + 1, j, labels(i + 1, j)));
            end
            if (j < M)
                energy = energy + ...
                    get_pp(i + N * (j - 1), i + N * j, ...
                    labels(i, j), labels(i, j + 1), intensity);
%                     abs(intensity(i, j, labels(i, j)) - intensity(i, j, labels(i, j + 1))) + ...
%                     abs(intensity(i, j + 1, labels(i, j)) - intensity(i, j + 1, labels(i, j + 1)));
            end
        end
    end
end