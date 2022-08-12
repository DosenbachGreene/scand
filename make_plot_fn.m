function f = make_plot_fn(sat, gpa, school)
% Closure used in example.m to generate scatter plots.

% Fixed effect estimates of beta (without clustering).
Bfixed = pinv([ones(size(sat,1),1), gpa])*sat;

% Closure that stores the data for sat and gpa, and generates plots based
% on fixed effect parameters B and random effect parameters u.
function plot_clusters(B, u)
    xrange = [0.9, 4.1];
    if nargin == 1
        figure; hold on;
        scatter_h = scatter(gpa, sat, 'k');
        xlabel('Grade Point Average'); ylabel('SAT Score');
        xlim(xrange); ylim([350,1650]);
        line(xrange, B(1) + xrange.*B(2), 'Color', scatter_h.CData, 'LineStyle', '--');
    else
        figure; hold on;
        scatter_h1 = scatter(gpa(logical(school(:,1))), sat(logical(school(:,1))));
        scatter_h2 = scatter(gpa(logical(school(:,2))), sat(logical(school(:,2))));
        scatter_h3 = scatter(gpa(logical(school(:,3))), sat(logical(school(:,3))));
        xlabel('Grade Point Average'); ylabel('SAT Score');
        xlim(xrange); ylim([350,1650]);
        line_fixed = line(xrange, Bfixed(1) + xrange.*Bfixed(2), 'Color', 'black', 'LineStyle', '--');
        line_marginal = line(xrange, B(1) + xrange.*B(2), 'Color', 'black');
        line(xrange, B(1) + u(1) + xrange.*(B(2) + u(2)), 'Color', scatter_h1.CData);
        line(xrange, B(1) + u(3) + xrange.*(B(2) + u(4)), 'Color', scatter_h2.CData);
        line(xrange, B(1) + u(5) + xrange.*(B(2) + u(6)), 'Color', scatter_h3.CData);
        legend([scatter_h3, scatter_h2, scatter_h1, line_marginal, line_fixed], 'Burroughs', 'Gateway', 'Jennings', 'Marginal', 'Fixed Only', 'Location', 'Northwest');
    end

end

% Return the closure.
f = @plot_clusters;

end
