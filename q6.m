n = 51
mus = linspace(-20, 30, n);
vals = bvp(mus(1));
for i = 2:n
    vals = vertcat(vals, bvp(mus(i)));
end

hold on

vals(:,2)
vals(:,3)

%u(1/4)
plot(mus, vals(:,1));
text(mus(1) + 2, vals(1,1) + 0.1, 'u(1/4)');
%u(1/2)
plot(mus, vals(:,2));
text(mus(1) + 2, vals(1,2) + 0.2, 'u(1/2)');
%u(3/4)
plot(mus, vals(:,3));
text(mus(1) + 2, vals(1,3) + 0.3, 'u(3/4)');

hold off
