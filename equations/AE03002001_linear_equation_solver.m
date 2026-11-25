% linear_equation_solver.m
% ============================================
% Interactive equation solver in GNU Octave
% The program asks for the number of unknowns and the degree of the
% equation, then:
%   - degree 1 (linear): solves a system of linear equations A*x = b
%     using several methods, classifies it (unique / infinite / no
%     solution) and draws many kinds of plots.
%   - degree >= 2 with ONE unknown: solves the polynomial equation
%     a_d*x^d + ... + a_1*x + a_0 = 0 and plots it.
%
% Run with:  octave linear_equation_solver.m
% Tip: if the figure windows do not respond while the program
% waits for input, run it from the Octave prompt instead:
%   >> linear_equation_solver
% ============================================
1;  % Marks this file as a script (functions are defined below)

% ============================================
% Input helpers
% ============================================

% --- Reads a positive integer ---
function n = read_positive_int(prompt_text)
  while true
    n = str2double(strtrim(input(prompt_text, 's')));
    if isscalar(n) && isfinite(n) && n >= 1 && n == fix(n)
      return;
    end
    disp('Invalid input. Please enter a positive integer.');
  end
end

% --- Reads a finite real number (coefficient / menu option) ---
function v = read_number(prompt_text)
  while true
    v = str2double(strtrim(input(prompt_text, 's')));
    if isscalar(v) && isfinite(v)
      return;
    end
    disp('Invalid input. Please enter a number (decimals allowed, e.g. -2.5).');
  end
end

% --- Asks whether to save a figure as a PNG file ---
function offer_save(fig, filename)
  answer = lower(strtrim(input(sprintf('Save figure as %s? (y/n): ', filename), 's')));
  if strcmp(answer, 'y')
    print(fig, filename, '-dpng', '-r150');
    printf('Saved: %s\n', filename);
  end
end

% ============================================
% Linear systems: input, printing and analysis
% ============================================

% --- Reads the coefficients of m equations with n unknowns ---
function [A, b] = read_system(m, n)
  printf('\nEnter the coefficients of each equation: a1*x1 + ... + an*xn = b\n');
  A = zeros(m, n);
  b = zeros(m, 1);
  for i = 1:m
    printf('\n--- Equation %d ---\n', i);
    for j = 1:n
      A(i,j) = read_number(sprintf('  Coefficient of x%d: ', j));
    end
    b(i) = read_number('  Constant term (right-hand side): ');
  end
  printf('\nSystem entered (%d equations, %d unknowns):\n', m, n);
  print_system(A, b);
end

% --- Prints the system in a readable form, e.g.  2*x1 - 3*x2 = 5 ---
function print_system(A, b)
  [m, n] = size(A);
  for i = 1:m
    s = '';
    for j = 1:n
      a = A(i,j);
      if a == 0
        continue;
      end
      if isempty(s)
        s = sprintf('%g*x%d', a, j);
      else
        if a < 0
          sgn = '-';
        else
          sgn = '+';
        end
        s = [s, sprintf(' %s %g*x%d', sgn, abs(a), j)];
      end
    end
    if isempty(s)
      s = '0';
    end
    printf('  %s = %g\n', s, b(i));
  end
end

% --- Classifies the system with the Rouche-Capelli theorem ---
% status is 'unique', 'infinite' or 'inconsistent'
function [status, r, ra] = classify_system(A, b)
  r  = rank(A);
  ra = rank([A, b]);
  if r < ra
    status = 'inconsistent';
  elseif r == columns(A)
    status = 'unique';
  else
    status = 'infinite';
  end
end

% --- Prints rank, determinant and the type of solution ---
function analyze_system(A, b)
  [m, n] = size(A);
  [status, r, ra] = classify_system(A, b);
  printf('Equations: %d, unknowns: %d\n', m, n);
  printf('rank(A) = %d, rank([A|b]) = %d\n', r, ra);
  if m == n
    printf('det(A) = %.6g\n', det(A));
    printf('cond(A) = %.6g\n', cond(A));
  end
  switch status
    case 'unique'
      disp('Rouche-Capelli: rank(A) = rank([A|b]) = number of unknowns.');
      disp('=> The system has exactly ONE solution.');
    case 'infinite'
      disp('Rouche-Capelli: rank(A) = rank([A|b]) < number of unknowns.');
      printf('=> The system has INFINITELY many solutions (%d free variable(s)).\n', n - r);
    case 'inconsistent'
      disp('Rouche-Capelli: rank(A) < rank([A|b]).');
      disp('=> The system has NO solution (a least-squares solution can still be computed).');
  end
end

% --- Shows A*x, b and the residual of a computed solution ---
function verify_solution(A, b, x)
  if isempty(x)
    disp('No solution available yet. Use the "Solve" option first.');
    return;
  end
  res = A * x - b;
  disp('Verification (A*x compared with b):');
  for i = 1:rows(A)
    printf('  Eq %d: A*x = %.6g, b = %.6g, residual = %.3g\n', i, A(i,:) * x, b(i), res(i));
  end
  printf('Residual norm ||A*x - b|| = %.3g\n', norm(res));
end

% ============================================
% Linear systems: solution methods
% ============================================

% --- Asks for a method and solves the system. Returns [] if it fails ---
function x = solve_system(A, b)
  x = [];
  [m, n] = size(A);
  [status, r, ~] = classify_system(A, b);
  invertible = (m == n) && (r == n);   % square with non-zero determinant

  disp(' ');
  disp('--- Solution methods ---');
  disp(' 1) Backslash operator (A \ b)');
  disp(' 2) Gauss-Jordan elimination (rref)');
  disp(' 3) Cramer''s rule (square system, det(A) ~= 0)');
  disp(' 4) Inverse matrix: x = inv(A) * b (square system, det(A) ~= 0)');
  disp(' 5) LU decomposition (square system, det(A) ~= 0)');
  disp(' 6) Pseudoinverse: x = pinv(A) * b (least squares / minimum norm)');
  k = read_number('Method: ');

  label = 'Solution';
  switch k
    case 1
      if ~strcmp(status, 'unique')
        printf('Not applicable: this method needs a unique solution (status: %s). Try method 2 or 6.\n', status);
        return;
      end
      x = A \ b;

    case 2
      Ab = [A, b];
      [R, piv] = rref(Ab);
      disp('Reduced row echelon form of [A|b]:');
      disp(R);
      if strcmp(status, 'inconsistent')
        disp('The system has no solution. Use method 6 for a least-squares solution.');
        return;
      end
      x = zeros(n, 1);
      for q = 1:numel(piv)
        x(piv(q)) = R(q, end);   % free variables are set to 0
      end
      free = setdiff(1:n, piv);
      if ~isempty(free)
        label = 'Particular solution (free variables = 0)';
        printf('Free variable(s): ');
        printf('x%d ', free);
        printf('\n');
        for t = 1:numel(free)
          v = zeros(n, 1);
          v(free(t)) = 1;
          for q = 1:numel(piv)
            v(piv(q)) = -R(q, free(t));
          end
          printf('Direction vector v%d (for free variable x%d):\n', t, free(t));
          disp(v);
        end
        disp('General solution: x = xp + t1*v1 + t2*v2 + ...   (t1, t2, ... any real numbers)');
      end

    case 3
      if ~invertible
        disp('Not applicable: the system must be square with det(A) ~= 0.');
        return;
      end
      dA = det(A);
      printf('det(A) = %.6g\n', dA);
      x = zeros(n, 1);
      for i = 1:n
        Ai = A;
        Ai(:, i) = b;   % replace column i by the constant terms
        printf('det(A_%d) = %.6g\n', i, det(Ai));
        x(i) = det(Ai) / dA;
      end

    case 4
      if ~invertible
        disp('Not applicable: the system must be square with det(A) ~= 0.');
        return;
      end
      Ainv = inv(A);
      disp('inv(A) ='); disp(Ainv);
      x = Ainv * b;

    case 5
      if ~invertible
        disp('Not applicable: the system must be square with det(A) ~= 0.');
        return;
      end
      [L, U, P] = lu(A);   % P*A = L*U
      disp('L ='); disp(L);
      disp('U ='); disp(U);
      disp('P ='); disp(P);
      y = L \ (P * b);     % forward substitution
      x = U \ y;           % backward substitution

    case 6
      x = pinv(A) * b;
      if strcmp(status, 'inconsistent')
        label = 'Least-squares solution';
      elseif strcmp(status, 'infinite')
        label = 'Minimum-norm solution';
      end

    otherwise
      disp('Invalid method.');
      return;
  end

  printf('\n%s:\n', label);
  for i = 1:n
    printf('  x%d = %.6g\n', i, x(i));
  end
  printf('Residual norm ||A*x - b|| = %.3g\n', norm(A * x - b));
end

% ============================================
% Linear systems: plots
% ============================================

% --- Geometric view: lines (1-2 unknowns) or planes (3 unknowns) ---
function draw_geometry(A, b, x)
  [m, n] = size(A);
  x = x(:);
  w = max(3, 1.5 * max(abs(x)));   % half-width of the view
  cols = lines(m);
  hs = [];
  labels = {};
  hold on;

  if n == 1
    xs = linspace(x(1) - w, x(1) + w, 200);
    for i = 1:m
      h = plot(xs, A(i,1) * xs - b(i), 'Color', cols(i,:), 'LineWidth', 1.5);
      hs(end+1) = h;
      labels{end+1} = sprintf('Eq %d', i);
    end
    plot(xs([1 end]), [0 0], 'k--');
    hp = plot(x(1), 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    xlabel('x1'); ylabel('a*x1 - b');
    title('Each equation as a*x1 - b (zero crossing = solution)');
    grid on;

  elseif n == 2
    xs = linspace(x(1) - w, x(1) + w, 200);
    for i = 1:m
      a1 = A(i,1); a2 = A(i,2);
      if a2 ~= 0
        h = plot(xs, (b(i) - a1 * xs) / a2, 'Color', cols(i,:), 'LineWidth', 1.5);
      elseif a1 ~= 0
        h = plot(b(i) / a1 * [1 1], [x(2) - w, x(2) + w], 'Color', cols(i,:), 'LineWidth', 1.5);
      else
        continue;
      end
      hs(end+1) = h;
      labels{end+1} = sprintf('Eq %d', i);
    end
    hp = plot(x(1), x(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    axis([x(1) - w, x(1) + w, x(2) - w, x(2) + w]);
    xlabel('x1'); ylabel('x2');
    title('Lines of the equations');
    grid on;

  else   % n == 3
    rng = {linspace(x(1) - w, x(1) + w, 15), ...
           linspace(x(2) - w, x(2) + w, 15), ...
           linspace(x(3) - w, x(3) + w, 15)};
    for i = 1:m
      a = A(i,:);
      k = find(a ~= 0, 1, 'last');   % unknown that is solved from the plane equation
      if isempty(k)
        continue;
      end
      others = setdiff(1:3, k);
      [P, Q] = meshgrid(rng{others(1)}, rng{others(2)});
      S = (b(i) - a(others(1)) * P - a(others(2)) * Q) / a(k);
      coords = cell(1, 3);
      coords{others(1)} = P;
      coords{others(2)} = Q;
      coords{k} = S;
      h = surf(coords{1}, coords{2}, coords{3}, 'FaceColor', cols(i,:), 'EdgeColor', 'none');
      try
        set(h, 'FaceAlpha', 0.5);
      catch
      end
      hs(end+1) = h;
      labels{end+1} = sprintf('Eq %d', i);
    end
    hp = plot3(x(1), x(2), x(3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    axis([x(1) - w, x(1) + w, x(2) - w, x(2) + w, x(3) - w, x(3) + w]);
    xlabel('x1'); ylabel('x2'); zlabel('x3');
    title('Planes of the equations');
    view(-37.5, 30);
    rotate3d on;
    grid on;
  end

  hold off;
  hs(end+1) = hp;
  labels{end+1} = 'Computed x';
  try
    legend(hs, labels);
  catch
  end
end

% --- Bar chart of the solution vector ---
function draw_solution_bar(x)
  x = x(:);
  n = numel(x);
  bar(1:n, x);
  title('Solution vector'); xlabel('Unknown'); ylabel('Value');
  grid on;
  if n <= 20
    set(gca, 'XTick', 1:n, 'XTickLabel', arrayfun(@(k) sprintf('x%d', k), 1:n, 'UniformOutput', false));
  end
  if n <= 12
    for k = 1:n
      if x(k) >= 0
        va = 'bottom';
      else
        va = 'top';
      end
      text(k, x(k), sprintf('%.3g', x(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', va);
    end
  end
end

% --- Stem chart of the solution vector ---
function draw_solution_stem(x)
  x = x(:);
  n = numel(x);
  stem(1:n, x, 'filled');
  title('Solution vector (stem)'); xlabel('Unknown'); ylabel('Value');
  grid on;
  if n <= 20
    set(gca, 'XTick', 1:n, 'XTickLabel', arrayfun(@(k) sprintf('x%d', k), 1:n, 'UniformOutput', false));
  end
end

% --- Grouped bars: A*x versus b for every equation ---
function draw_ax_vs_b(A, b, x)
  bar([A * x(:), b(:)]);
  legend({'A*x', 'b'});
  title('A*x compared with b'); xlabel('Equation'); ylabel('Value');
  grid on;
end

% --- Residual A*x - b for every equation ---
function draw_residual(A, b, x)
  res = A * x(:) - b(:);
  bar(res);
  title(sprintf('Residual per equation (norm = %.3g)', norm(res)));
  xlabel('Equation'); ylabel('A*x - b');
  grid on;
end

% --- Heat map of the augmented matrix [A|b] ---
function draw_heatmap(A, b)
  [m, n] = size(A);
  M = [A, b(:)];
  imagesc(M);
  colorbar;
  hold on;
  plot([n + 0.5, n + 0.5], [0.5, m + 0.5], 'k-', 'LineWidth', 2);   % separates A from b
  hold off;
  title('Augmented matrix [A | b]'); xlabel('Column'); ylabel('Equation');
  if m <= 10 && n + 1 <= 10
    for i = 1:m
      for j = 1:n + 1
        text(j, i, sprintf('%.3g', M(i,j)), 'HorizontalAlignment', 'center', ...
             'Color', 'w', 'FontWeight', 'bold');
      end
    end
  end
  if n + 1 <= 20
    labels = [arrayfun(@(k) sprintf('x%d', k), 1:n, 'UniformOutput', false), {'b'}];
    set(gca, 'XTick', 1:n + 1, 'XTickLabel', labels);
  end
  if m <= 20
    set(gca, 'YTick', 1:m);
  end
end

% --- Smooth surface of the augmented matrix [A|b] ---
% The matrix is interpolated onto a fine grid; the original values are black dots.
function draw_aug_surface(A, b)
  M = [A, b(:)];
  [r, c] = size(M);
  [X, Y] = meshgrid(1:c, 1:r);
  [Xi, Yi] = meshgrid(linspace(1, c, 60), linspace(1, r, 60));
  if min(r, c) >= 3
    method = 'spline';
  else
    method = 'linear';
  end
  Zi = interp2(X, Y, M, Xi, Yi, method);
  surf(Xi, Yi, Zi);
  shading interp;
  hold on;
  plot3(X(:), Y(:), M(:), 'ko', 'MarkerFaceColor', 'k');
  hold off;
  colorbar;
  title('Surface of [A | b]'); xlabel('Column'); ylabel('Equation'); zlabel('Value');
  view(-37.5, 30);
  rotate3d on;
end

% --- 3D bars of the augmented matrix [A|b] (stem3 if bar3 is missing) ---
function draw_aug_bars3d(A, b)
  M = [A, b(:)];
  [r, c] = size(M);
  if exist('bar3') > 0
    bar3(M);
  else
    disp('Note: bar3 is not available in this Octave version; using stem3 instead.');
    [X, Y] = meshgrid(1:c, 1:r);
    stem3(X(:), Y(:), M(:), 'filled');
  end
  title('3D bars of [A | b]'); xlabel('Column'); ylabel('Equation'); zlabel('Value');
  view(-37.5, 30);
  rotate3d on;
end

% --- Singular values of A ---
function draw_singular(A)
  s = svd(A);
  bar(s);
  title(sprintf('Singular values of A (rank = %d)', rank(A)));
  xlabel('Index'); ylabel('Singular value');
  grid on;
end

% --- Menu of plots for a linear system ---
function plot_linear(A, b, x)
  [m, n] = size(A);
  if isempty(x)
    x = pinv(A) * b;
    disp('Note: no solution computed yet; the plots use pinv(A)*b (least squares / minimum norm).');
  end

  names = {'geometry', 'solution_bar', 'solution_stem', 'ax_vs_b', 'residual', ...
           'heatmap', 'surface', 'bars3d', 'singular_values', 'dashboard'};
  disp(' ');
  disp('--- Plot types ---');
  disp(' 1) Geometric view (lines / planes, only for 1, 2 or 3 unknowns)');
  disp(' 2) Solution vector (bar)');
  disp(' 3) Solution vector (stem)');
  disp(' 4) A*x compared with b (grouped bars)');
  disp(' 5) Residual per equation (bar)');
  disp(' 6) Augmented matrix [A|b] (heat map)');
  disp(' 7) Augmented matrix [A|b] (smooth surface, needs 2+ equations)');
  disp(' 8) Augmented matrix [A|b] (3D bars)');
  disp(' 9) Singular values of A (bar)');
  disp('10) Dashboard (several plots in one figure)');
  k = read_number('Plot type: ');

  if ~(k >= 1 && k <= numel(names) && k == fix(k))
    disp('Invalid plot type.');
    return;
  end
  if k == 1 && n > 3
    disp('The geometric view is only available for 1, 2 or 3 unknowns.');
    return;
  end
  if k == 7 && m < 2
    disp('Error: the surface needs at least 2 equations.');
    return;
  end

  if k == 10
    fig = figure('Name', 'Linear system dashboard', 'Position', [100 100 1300 750]);
  else
    fig = figure('Name', names{k}, 'Position', [100 100 800 600]);
  end

  switch k
    case 1
      draw_geometry(A, b, x);
    case 2
      draw_solution_bar(x);
    case 3
      draw_solution_stem(x);
    case 4
      draw_ax_vs_b(A, b, x);
    case 5
      draw_residual(A, b, x);
    case 6
      draw_heatmap(A, b);
    case 7
      draw_aug_surface(A, b);
    case 8
      draw_aug_bars3d(A, b);
    case 9
      draw_singular(A);
    case 10
      subplot(2, 3, 1);
      if n <= 3
        draw_geometry(A, b, x);
      else
        draw_solution_stem(x);
      end
      subplot(2, 3, 2); draw_solution_bar(x);
      subplot(2, 3, 3); draw_ax_vs_b(A, b, x);
      subplot(2, 3, 4); draw_residual(A, b, x);
      subplot(2, 3, 5); draw_heatmap(A, b);
      subplot(2, 3, 6); draw_singular(A);
  end
  drawnow;
  offer_save(fig, sprintf('linear_system_%s.png', names{k}));
end

% ============================================
% Linear system mode (degree 1)
% ============================================

% --- Menu loop for a linear system. Returns true to enter a new problem ---
function again = linear_mode(n)
  m = read_positive_int('Number of equations: ');
  [A, b] = read_system(m, n);
  x = [];

  while true
    disp(' ');
    disp('========== MENU (linear system) ==========');
    disp(' 1) Show the system');
    disp(' 2) Analyze the system (rank, determinant, type of solution)');
    disp(' 3) Solve the system (choose a method)');
    disp(' 4) Verify the solution (residual)');
    disp(' 5) Plot (choose the plot type)');
    disp(' 6) Enter a new problem');
    disp(' 0) Exit');
    option = read_number('Choose an option: ');

    switch option
      case 1
        print_system(A, b);

      case 2
        analyze_system(A, b);

      case 3
        new_x = solve_system(A, b);
        if ~isempty(new_x)
          x = new_x;
        end

      case 4
        verify_solution(A, b, x);

      case 5
        plot_linear(A, b, x);

      case 6
        again = true;
        return;

      case 0
        again = false;
        return;

      otherwise
        disp('Invalid option. Please choose a number from the menu.');
    end
  end
end

% ============================================
% Polynomial equations (degree >= 2, one unknown)
% ============================================

% --- Text form of a polynomial equation, e.g.  2*x^2 - 3*x + 1 = 0 ---
function s = poly_to_string(p)
  d = numel(p) - 1;
  s = '';
  for k = 1:numel(p)
    c = p(k);
    e = d - k + 1;   % exponent of this term
    if c == 0
      continue;
    end
    if e == 0
      term = sprintf('%g', abs(c));
    elseif e == 1
      term = sprintf('%g*x', abs(c));
    else
      term = sprintf('%g*x^%d', abs(c), e);
    end
    if isempty(s)
      if c < 0
        s = ['-', term];
      else
        s = term;
      end
    else
      if c < 0
        s = [s, ' - ', term];
      else
        s = [s, ' + ', term];
      end
    end
  end
  if isempty(s)
    s = '0';
  end
  s = [s, ' = 0'];
end

% --- Prints the roots of the polynomial and checks them ---
function solve_polynomial(p)
  d = numel(p) - 1;
  r = roots(p);
  if d == 2
    D = p(2)^2 - 4 * p(1) * p(3);
    printf('Discriminant = %.6g  ->  ', D);
    if D > 0
      disp('two distinct real roots.');
    elseif D == 0
      disp('one double real root.');
    else
      disp('two complex conjugate roots.');
    end
  end
  printf('Roots (%d, counted with multiplicity):\n', d);
  for k = 1:numel(r)
    if abs(imag(r(k))) < 1e-6 * max(1, abs(r(k)))
      printf('  x%d = %.6g\n', k, real(r(k)));
    else
      if imag(r(k)) < 0
        sgn = '-';
      else
        sgn = '+';
      end
      printf('  x%d = %.6g %s %.6gi\n', k, real(r(k)), sgn, abs(imag(r(k))));
    end
  end
  printf('Verification: max |f(root)| = %.3g\n', max(abs(polyval(p, r))));
end

% --- Real parts of the (numerically) real roots ---
function rr = real_roots(r)
  tol = 1e-6 * max(1, abs(r));
  rr = real(r(abs(imag(r)) < tol));
end

% --- Vertical limit that shows the interesting part of a polynomial ---
function ymax = poly_ymax(p, R)
  ymax = 1;
  if numel(p) >= 3
    rc = real_roots(roots(polyder(p)));   % critical points
    rc = rc(abs(rc) <= 1.5 * R);
    if ~isempty(rc)
      ymax = max(1, max(abs(polyval(p, rc))));
    end
  end
  ymax = 1.3 * ymax;
end

% --- Curve of the polynomial with its real roots ---
function draw_poly_curve(p, r)
  R = max(1, max(abs(r)));
  xs = linspace(-1.5 * R, 1.5 * R, 600);
  plot(xs, polyval(p, xs), 'b-', 'LineWidth', 1.5);
  hold on;
  plot(xs([1 end]), [0 0], 'k-');
  rr = real_roots(r);
  plot(rr, zeros(size(rr)), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
  hold off;
  ylim([-1, 1] * poly_ymax(p, R));
  grid on;
  xlabel('x'); ylabel('f(x)');
  title(['f(x): ', poly_to_string(p)]);
end

% --- Roots in the complex plane ---
function draw_roots_plane(r)
  R = 1.3 * max(1, max(abs(r)));
  plot(real(r), imag(r), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
  hold on;
  plot([-R R], [0 0], 'k-');
  plot([0 0], [-R R], 'k-');
  hold off;
  axis([-R R -R R]);
  axis square;
  grid on;
  xlabel('Real part'); ylabel('Imaginary part');
  title('Roots in the complex plane');
end

% --- f(x), f'(x) and f''(x) ---
function draw_poly_derivatives(p)
  r = roots(p);
  R = max(1, max(abs(r)));
  xs = linspace(-1.5 * R, 1.5 * R, 600);
  dp = polyder(p);
  ddp = polyder(dp);
  h = plot(xs, polyval(p, xs), 'b-', xs, polyval(dp, xs), 'r--', xs, polyval(ddp, xs), 'g:');
  set(h, 'LineWidth', 1.5);
  legend(h, {'f(x)', 'f''(x)', 'f''''(x)'});
  ymax = max([poly_ymax(p, R), poly_ymax(dp, R), poly_ymax(ddp, R)]);
  ylim([-ymax, ymax]);
  grid on;
  xlabel('x'); ylabel('Value');
  title('Polynomial and its first two derivatives');
end

% --- Menu of plots for a polynomial equation ---
function plot_polynomial(p)
  r = roots(p);
  names = {'curve', 'complex_plane', 'curve_and_complex_plane', 'derivatives'};
  disp(' ');
  disp('--- Plot types ---');
  disp(' 1) Curve f(x) with its real roots');
  disp(' 2) Roots in the complex plane');
  disp(' 3) Both plots side by side');
  disp(' 4) f(x) and its first two derivatives');
  k = read_number('Plot type: ');

  if ~(k >= 1 && k <= numel(names) && k == fix(k))
    disp('Invalid plot type.');
    return;
  end

  if k == 3
    fig = figure('Name', 'Polynomial', 'Position', [100 100 1200 500]);
  else
    fig = figure('Name', names{k}, 'Position', [100 100 800 600]);
  end

  switch k
    case 1
      draw_poly_curve(p, r);
    case 2
      draw_roots_plane(r);
    case 3
      subplot(1, 2, 1); draw_poly_curve(p, r);
      subplot(1, 2, 2); draw_roots_plane(r);
    case 4
      draw_poly_derivatives(p);
  end
  drawnow;
  offer_save(fig, sprintf('polynomial_%s.png', names{k}));
end

% --- Menu loop for a polynomial equation. Returns true to enter a new problem ---
function again = polynomial_mode(d)
  printf('\nEnter the coefficients of the equation of degree %d (from highest to lowest power).\n', d);
  p = zeros(1, d + 1);
  for k = 1:d + 1
    e = d - k + 1;
    while true
      p(k) = read_number(sprintf('  Coefficient of x^%d: ', e));
      if k == 1 && p(k) == 0
        printf('The leading coefficient cannot be zero in an equation of degree %d.\n', d);
      else
        break;
      end
    end
  end
  printf('\nEquation entered:\n  %s\n', poly_to_string(p));

  while true
    disp(' ');
    disp('========== MENU (polynomial equation) ==========');
    disp(' 1) Show the equation');
    disp(' 2) Solve the equation (roots)');
    disp(' 3) Plot (choose the plot type)');
    disp(' 4) Enter a new problem');
    disp(' 0) Exit');
    option = read_number('Choose an option: ');

    switch option
      case 1
        printf('  %s\n', poly_to_string(p));

      case 2
        solve_polynomial(p);

      case 3
        plot_polynomial(p);

      case 4
        again = true;
        return;

      case 0
        again = false;
        return;

      otherwise
        disp('Invalid option. Please choose a number from the menu.');
    end
  end
end

% ============================================
% Main program
% ============================================
disp('=== Equation Solver ===');
again = true;
while again
  disp(' ');
  n = read_positive_int('Number of unknowns: ');
  d = read_positive_int('Degree of the equation (1 = linear): ');

  if d == 1
    again = linear_mode(n);
  elseif n == 1
    again = polynomial_mode(d);
  else
    disp('Equations of degree 2 or higher are only supported with ONE unknown.');
    disp('Use degree 1 for several unknowns, or 1 unknown for a higher degree.');
  end
end
disp('Goodbye!');
