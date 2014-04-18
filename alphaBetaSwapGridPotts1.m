function [labels, energy, time] = alphaBetaSwapGridPotts1(unary, vertC, horC, metric, labels)
    success = 1;
    N = size( unary,1);
    M = size( unary, 2);
    K = size( unary, 3);
%     load('labels.mat');
    %labels = unidrnd(K, [ N M]);
    unaryR = zeros(N,M * K);
    for i = 1 : K
        unaryR(:, 1 + (i - 1) * M : i * M) = unary(:, :, i);
    end
   while(success)
        success = 0;
        for alpha = 1 : K - 1
            for beta  = alpha + 1 : K
                P_ab = (labels ==alpha) | (labels == beta);
                t = unary(:,:,alpha);
                s = unary(:, :, beta);
                for i = 1 : N
                    for j = 1 : M
                        if P_ab(i,j) == 1
                            if( i  > 1)
                                t(i,j) = t(i, j) + vertC( i - 1, j ) *  ( ~ P_ab(i - 1, j )) *  metric(labels(i - 1, j) +  K *  (alpha - 1)); 
                                s(i, j) = s(i, j) + vertC( i - 1, j) *  ( ~ P_ab(i - 1, j )) *  metric(labels(i - 1, j) +  K *  (beta - 1));
                                %t(i,j) = t(i,j) + vertC(i - 1, j) *( ( ~ P_ab(i - 1, j )) *  metric(alpha +  K *  (labels(i - 1, j ) - 1)) - P_ab(i - 1, j ) *(metric(beta +  K *  ( beta - 1))  - metric(alpha +  K *  ( alpha - 1))) );
                                %s(i,j) = s(i,j) +  vertC(i - 1, j ) * (( ~ P_ab(i - 1, j )) *  metric(beta + K  *  (labels(i - 1, j) - 1))) ;
                            end
                            if( i <  N )
                                t(i,j) = t(i,j) + vertC(i, j)  * ( ( ~ P_ab(i + 1, j )) *  metric(alpha +  K * ( labels(i + 1, j ) - 1)) +  P_ab(i + 1, j ) * metric(alpha +  K *  ( alpha - 1)) );
                                s(i,j) = s(i,j) +  vertC(i , j) * ( ( ~P_ab(i + 1, j))  *  metric(beta +  K *  (labels(i + 1, j ) - 1))  + P_ab(i + 1, j) *  metric(beta +  K *  ( beta - 1))) ;
                            end
                            if( j > 1 )
                                t(i, j) = t(i, j) + horC( i, j - 1) * ( ~P_ab(i , j - 1 ))  * metric(alpha + K* (labels(i, j - 1) - 1)) ; 
                                s(i,j) = s(i,j) + ( ~P_ab(i , j - 1 )) * horC(i , j - 1) *  metric(beta +  K * ( labels(i , j - 1) - 1)); 
                                %t(i,j) = t(i,j) + horC(i, j - 1) *(( ~P_ab(i , j - 1 ))  * metric(alpha + K* (labels(i, j - 1) - 1)) - P_ab(i , j - 1 ) * (metric(beta +  K *  ( beta - 1))) - metric(alpha +  K *  ( alpha - 1)));
                                %s(i,j) = s(i,j) + ( ~P_ab(i , j - 1 )) * horC(i , j - 1) *  metric(beta +  K * ( labels(i , j - 1) - 1));
                            end
                            if( j < M)
                                t(i,j) = t(i,j) +  horC(i, j )  * (( ~P_ab(i , j + 1 ))  * metric(alpha + K *  (labels(i, j + 1) - 1))  +  P_ab(i , j + 1 ) * (metric(alpha +  K *  (alpha - 1)) ));
                                s(i,j) = s(i,j) +  horC(i , j )  * (( ~P_ab(i , j + 1 )) *  metric(beta +  K * ( labels(i , j + 1) - 1)) + P_ab(i , j + 1 ) * (metric(beta +  K *  (beta- 1)) ));
                            end
                                
                        end
                    end    
                end
                
                %binary potencial
                
                vertex_hor = P_ab(:, 1 : end - 1) .* P_ab(: , 2 : end);
                vertex_vert = P_ab(1 : end - 1, :) .* P_ab(2 : end, :);
                theta_ab_hor = horC * (metric(alpha + K * (beta - 1)) - metric(alpha+ K * (alpha - 1)));
                theta_ab_vert = vertC * (metric(alpha + K * (beta - 1)) - metric(alpha+ K * (alpha - 1)));
                
                theta_ab_hor = theta_ab_hor(vertex_hor > 0);
                theta_ab_hor = theta_ab_hor(:);
                
                theta_ab_vert = theta_ab_vert(vertex_vert > 0);
                theta_ab_vert = theta_ab_vert(:);
                
                theta_ba_hor = horC * (metric(beta + K  *  (alpha - 1)) -    metric(beta+ K * (beta - 1)));
                theta_ba_vert = vertC * (metric(beta + K  *  (alpha - 1)) -    metric(beta+ K * (beta - 1)));
                
                theta_ba_hor = theta_ba_hor(vertex_hor > 0);
                theta_ba_hor = theta_ba_hor(:);
                
                theta_ba_vert = theta_ba_vert(vertex_vert > 0);
                theta_ba_vert = theta_ba_vert(:);
                
                %test this moment
                t = t(P_ab);
                s = s(P_ab);
                terminal_weights = [t(:), s(:)];
                
                %vertex making
                f = sum(P_ab(:));
                P_ab1 = P_ab * 1;
                P_ab1(P_ab) = 1 : f;

                vertex_num = P_ab1;
                
                vertex_hor2 =  vertex_hor .* vertex_num(:, 2 : end );
                vertex_hor = vertex_hor .* vertex_num(:, 1 : end - 1);
                
                vertex_vert = vertex_vert .* vertex_num(1 : end - 1, :);
                vertex_hor2 = vertex_hor2(vertex_hor2 ~=  0);
                vertex_hor = vertex_hor(vertex_hor ~=  0);
                vertex_vert = vertex_vert(vertex_vert ~= 0);
                
                %edge weights
                
                p1 = [vertex_vert,  vertex_vert + 1, theta_ba_vert, theta_ab_vert];
                p2 = [vertex_hor, vertex_hor2, theta_ba_hor, theta_ab_hor];
                
                edge_weights=[p1; p2];
                
               	energy1 = energy_calculation( labels, unaryR, vertC, horC, metric, K);
                [cut, labelsN] = graphCutMex(terminal_weights,edge_weights);
                labelsN = (labelsN == 0) * beta + (labelsN ~= 0) * alpha;
                labels1 = labels;
                labels1(P_ab) = labelsN;
                energy2 = energy_calculation( labels1, unaryR, vertC, horC, metric, K);
                if (energy2 < energy1)
                    labels = labels1;
                    success = 1;
                end
           % (1 - x_i) * x_j
            end
        end
    end
end