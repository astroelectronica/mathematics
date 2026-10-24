% AE03001001_interactive_matrix_operations_plots.m
% ============================================
% astroelectronic@
% Claude Sonnet 5.5
% AE03001001
% ============================================
% Interactive matrix operations in GNU Octave
% The program asks for the size and the content of
% matrices A and B, then offers a menu of operations
% and lets the user visualize the matrices with many
% different plot types (surf, mesh, bar, contour, ...).
%
% Run with:  octave interactive_matrix_operations_plots.m
% Tip: if the figure windows do not respond while the program
% waits for input, run it from the Octave prompt instead:
%   >> interactive_matrix_operations_plots
% ============================================
1;  % Marks this file as a script (functions are defined below)

% --- Reads a positive integer (rows / columns) ---
function n = read_positive_int(prompt_text)
  while true
    n = str2double(strtrim(input(prompt_text, 's')));
    if isscalar(n) && isfinite(n) && n >= 1 && n == fix(n)
      return;
    end
    disp('Invalid input. Please enter a positive integer.');
  end
end

% --- Reads a finite real number (matrix element / menu option) ---
function v = read_number(prompt_text)
  while true
    v = str2double(strtrim(input(prompt_text, 's')));
    if isscalar(v) && isfinite(v)
      return;
    end
    disp('Invalid input. Please enter a number (decimals allowed, e.g. -2.5).');
  end
end

% --- Asks for the size and every element of a matrix ---
function M = read_matrix(name)
  printf('\n--- Matrix %s ---\n', name);
  rows = read_positive_int(sprintf('Number of rows of %s: ', name));
  cols = read_positive_int(sprintf('Number of columns of %s: ', name));
  M = zeros(rows, cols);
  for i = 1:rows
    for j = 1:cols
      M(i,j) = read_number(sprintf('  %s(%d,%d) = ', name, i, j));
    end
  end
  printf('\nMatrix %s (%d x %d):\n', name, rows, cols);
  disp(M);
end

% --- Lets the user choose between A and B ---
function [M, name] = pick_matrix(A, B)
  while true
    choice = upper(strtrim(input('Which matrix? (A/B): ', 's')));
    if strcmp(choice, 'A')
      M = A; name = 'A'; return;
    elseif strcmp(choice, 'B')
      M = B; name = 'B'; return;
    end
    disp('Invalid choice. Type A or B.');
  end
end

% --- Catalog of the available plot types ---
function [keys, descs] = plot_catalog()
  keys  = {'surf', 'surfl', 'surfc', 'mesh', 'meshc', 'waterfall', ...
           'contour', 'contour3', 'contourf', 'pcolor', 'imagesc', ...
           'bar', 'barh', 'bar3', 'ribbon', 'stem3', 'scatter3', ...
           'plot', 'stairs', 'area', 'quiver', 'hist'};
  descs = {'Solid surface', ...
           'Surface with lighting', ...
           'Surface + contour lines underneath', ...
           'Wireframe mesh', ...
           'Mesh + contour lines underneath', ...
           'Waterfall (one line per row)', ...
           'Contour lines (2D)', ...
           'Contour lines in 3D', ...
           'Filled contour lines', ...
           'Smooth color map (pcolor)', ...
           'Heat map (shows the matrix values)', ...
           'Grouped vertical bars (one series per column)', ...
           'Grouped horizontal bars', ...
           '3D bars', ...
           '3D ribbons (one per column)', ...
           '3D stems', ...
           '3D scatter colored by value', ...
           'Line plot (one line per column)', ...
           'Stairs plot (one series per column)', ...
           'Stacked area (one layer per column)', ...
           'Gradient vector field (arrows)', ...
           'Histogram of all the matrix values'};
end

% --- Function Octave must provide for each plot type (if not the same name) ---
% Returns the fallback type to use when the needed function is missing.
function ptype = resolve_plot_type(ptype)
  fallback = struct('bar3', 'stem3', 'surfl', 'surf', 'waterfall', 'mesh', ...
                    'ribbon', 'bar', 'contour3', 'contour', ...
                    'scatter3', 'stem3', 'area', 'bar');
  if isfield(fallback, ptype) && exist(ptype) == 0
    new_type = fallback.(ptype);
    printf('Note: "%s" is not available in this Octave version; using "%s" instead.\n', ...
           ptype, new_type);
    ptype = new_type;
  end
end

% --- True if the plot type needs at least 2 rows and 2 columns ---
function tf = needs_grid(ptype)
  tf = any(strcmp(ptype, {'surf','surfl','surfc','mesh','meshc','waterfall', ...
                          'contour','contour3','contourf','pcolor','quiver'}));
end

% --- True if the plot type is 3D (it can be rotated with the mouse) ---
function tf = is_3d(ptype)
  tf = any(strcmp(ptype, {'surf','surfl','surfc','mesh','meshc','waterfall', ...
                          'contour3','bar3','ribbon','stem3','scatter3'}));
end

% --- True if the plot type draws one series per column (2D charts) ---
function tf = is_series_chart(ptype)
  tf = any(strcmp(ptype, {'bar','barh','plot','stairs','area','hist'}));
end

% --- True if a colorbar makes sense for the plot type ---
function tf = uses_colorbar(ptype)
  tf = any(strcmp(ptype, {'surf','surfl','surfc','mesh','meshc','contour', ...
                          'contour3','contourf','pcolor','imagesc','scatter3'}));
end

% --- Shows the plot type menu and returns the chosen one ---
function ptype = choose_plot_type()
  [keys, descs] = plot_catalog();
  disp(' ');
  disp('--- Plot types ---');
  for k = 1:numel(keys)
    printf('%2d) %-10s %s\n', k, keys{k}, descs{k});
  end
  while true
    k = read_number('Plot type: ');
    if k >= 1 && k <= numel(keys) && k == fix(k)
      ptype = resolve_plot_type(keys{k});
      return;
    end
    disp('Invalid option. Choose a number from the list.');
  end
end

% --- Draws a matrix with the given plot type ---
% For surface / contour plots the matrix is interpolated onto a fine grid
% so the result is smooth. The original matrix values are shown as black dots.
function draw_plot(M, title_text, ptype)
  [r, c] = size(M);
  [X, Y] = meshgrid(1:c, 1:r);

  if needs_grid(ptype)
    if min(r, c) >= 3
      method = 'spline';   % smooth curves through the data points
    else
      method = 'linear';   % bilinear interpolation for very small matrices
    end
    n = 60;   % resolution of the smooth grid (n x n points)
    [Xi, Yi] = meshgrid(linspace(1, c, n), linspace(1, r, n));
    Zi = interp2(X, Y, M, Xi, Yi, method);
    % Waterfall: smooth only along the columns, keep one line per row
    [Xw, Yw] = meshgrid(linspace(1, c, n), 1:r);
    Zw = interp2(X, Y, M, Xw, Yw, method);
  end

  switch ptype
    case 'surf'
      surf(Xi, Yi, Zi); shading interp;
    case 'surfl'
      surfl(Xi, Yi, Zi); shading interp;
    case 'surfc'
      surfc(Xi, Yi, Zi); shading interp;
    case 'mesh'
      mesh(Xi, Yi, Zi);
    case 'meshc'
      meshc(Xi, Yi, Zi);
    case 'waterfall'
      waterfall(Xw, Yw, Zw);
    case 'contour'
      contour(Xi, Yi, Zi, 15);
    case 'contour3'
      contour3(Xi, Yi, Zi, 15);
    case 'contourf'
      contourf(Xi, Yi, Zi, 15);
    case 'pcolor'
      pcolor(Xi, Yi, Zi); shading interp;
    case 'imagesc'
      imagesc(M);
      if r <= 10 && c <= 10   % write the values inside the cells
        for i = 1:r
          for j = 1:c
            text(j, i, sprintf('%.3g', M(i,j)), 'HorizontalAlignment', 'center', ...
                 'Color', 'w', 'FontWeight', 'bold');
          end
        end
      end
    case 'bar'
      bar(M);
    case 'barh'
      barh(M);
    case 'bar3'
      bar3(M);
    case 'ribbon'
      ribbon(M);
    case 'stem3'
      stem3(X(:), Y(:), M(:), 'filled');
    case 'scatter3'
      scatter3(X(:), Y(:), M(:), 80, M(:), 'filled');
    case 'plot'
      plot(M, '-o', 'LineWidth', 1.5); grid on;
    case 'stairs'
      stairs(M, 'LineWidth', 1.5); grid on;
    case 'area'
      area(M);
    case 'quiver'
      [gx, gy] = gradient(M);   % gradient along columns (gx) and rows (gy)
      quiver(X, Y, gx, gy);
      axis tight;
    case 'hist'
      nbins = max(5, min(20, ceil(sqrt(numel(M)))));
      hist(M(:), nbins);
  end

  % Black dots with the original values on top of the smooth plots
  if any(strcmp(ptype, {'surf','surfl','surfc','mesh','meshc','contour3'}))
    hold on;
    plot3(X(:), Y(:), M(:), 'ko', 'MarkerFaceColor', 'k');
    hold off;
  elseif any(strcmp(ptype, {'contour','contourf','pcolor'}))
    hold on;
    plot(X(:), Y(:), 'ko', 'MarkerFaceColor', 'k');
    hold off;
  end

  if uses_colorbar(ptype)
    colorbar;
  end
  title(title_text);

  % Axis labels, legends and ticks
  if strcmp(ptype, 'hist')
    xlabel('Value'); ylabel('Frequency');
  elseif is_series_chart(ptype)
    if strcmp(ptype, 'barh')
      ylabel('Row'); xlabel('Value');
    else
      xlabel('Row'); ylabel('Value');
    end
    if c <= 8
      legend(arrayfun(@(k) sprintf('Col %d', k), 1:c, 'UniformOutput', false));
    end
  else
    xlabel('Column'); ylabel('Row');
    if is_3d(ptype)
      zlabel('Value');
    end
    if c <= 20 && r <= 20
      set(gca, 'XTick', 1:c, 'YTick', 1:r);
    end
  end

  if is_3d(ptype)
    view(-37.5, 30);   % 3D viewing angle
    rotate3d on;       % drag with the mouse to rotate
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

% --- Plot of one matrix ---
function plot_matrix(M, name)
  ptype = choose_plot_type();
  if needs_grid(ptype) && min(size(M)) < 2
    printf('Error: "%s" needs a matrix with at least 2 rows and 2 columns.\n', ptype);
    return;
  end
  fig = figure('Name', sprintf('%s of %s', ptype, name), 'Position', [100 100 800 600]);
  draw_plot(M, sprintf('%s of %s', ptype, name), ptype);
  drawnow;
  offer_save(fig, sprintf('matrix_%s_%s.png', name, ptype));
end

% --- Plots of A, B and the result of an operation ---
function plot_operation(A, B)
  disp('1) A + B   2) A - B   3) A * B   4) A .* B');
  op = read_number('Operation to plot: ');

  switch op
    case {1, 2, 4}
      if ~isequal(size(A), size(B))
        disp('Error: A and B must have the same size.');
        return;
      end
      if op == 1
        R = A + B; label = 'A + B';
      elseif op == 2
        R = A - B; label = 'A - B';
      else
        R = A .* B; label = 'A .* B';
      end
    case 3
      if columns(A) ~= rows(B)
        printf('Error: columns of A (%d) must equal rows of B (%d).\n', ...
               columns(A), rows(B));
        return;
      end
      R = A * B; label = 'A * B';
    otherwise
      disp('Invalid operation.');
      return;
  end

  ptype = choose_plot_type();
  if needs_grid(ptype) && min([size(A), size(B), size(R)]) < 2
    printf('Error: "%s" needs matrices with at least 2 rows and 2 columns.\n', ptype);
    return;
  end

  fig = figure('Name', sprintf('%s of %s', ptype, label), 'Position', [100 100 1300 450]);
  subplot(1, 3, 1); draw_plot(A, 'Matrix A', ptype);
  subplot(1, 3, 2); draw_plot(B, 'Matrix B', ptype);
  subplot(1, 3, 3); draw_plot(R, sprintf('Result: %s', label), ptype);
  drawnow;
  offer_save(fig, sprintf('matrix_operation_%s.png', ptype));
end

% --- Main program ---
disp('=== Interactive Matrix Operations ===');
A = read_matrix('A');
B = read_matrix('B');

running = true;
while running
  disp(' ');
  disp('========== MENU ==========');
  disp(' 1) Show matrices A and B');
  disp(' 2) A + B');
  disp(' 3) A - B');
  disp(' 4) A * B (matrix product)');
  disp(' 5) A .* B (element-wise product)');
  disp(' 6) Transpose');
  disp(' 7) Determinant');
  disp(' 8) Inverse');
  disp(' 9) Rank');
  disp('10) Trace');
  disp('11) Plot a matrix (choose the plot type)');
  disp('12) Plot A, B and an operation result (choose the plot type)');
  disp('13) Enter new matrices');
  disp(' 0) Exit');
  option = read_number('Choose an option: ');

  switch option
    case 1
      disp('Matrix A:'); disp(A);
      disp('Matrix B:'); disp(B);

    case 2
      if isequal(size(A), size(B))
        disp('A + B ='); disp(A + B);
      else
        disp('Error: A and B must have the same size.');
      end

    case 3
      if isequal(size(A), size(B))
        disp('A - B ='); disp(A - B);
      else
        disp('Error: A and B must have the same size.');
      end

    case 4
      if columns(A) == rows(B)
        disp('A * B ='); disp(A * B);
      else
        printf('Error: columns of A (%d) must equal rows of B (%d).\n', ...
               columns(A), rows(B));
      end

    case 5
      if isequal(size(A), size(B))
        disp('A .* B ='); disp(A .* B);
      else
        disp('Error: A and B must have the same size.');
      end

    case 6
      M = pick_matrix(A, B);
      disp('Transpose ='); disp(M');

    case 7
      M = pick_matrix(A, B);
      if rows(M) == columns(M)
        printf('Determinant = %.4f\n', det(M));
      else
        disp('Error: the matrix must be square.');
      end

    case 8
      M = pick_matrix(A, B);
      if rows(M) ~= columns(M)
        disp('Error: the matrix must be square.');
      elseif rank(M) < rows(M)
        disp('Error: the matrix is singular (no inverse exists).');
      else
        disp('Inverse ='); disp(inv(M));
      end

    case 9
      M = pick_matrix(A, B);
      printf('Rank = %d\n', rank(M));

    case 10
      M = pick_matrix(A, B);
      if rows(M) == columns(M)
        printf('Trace = %g\n', trace(M));
      else
        disp('Error: the matrix must be square.');
      end

    case 11
      [M, name] = pick_matrix(A, B);
      plot_matrix(M, name);

    case 12
      plot_operation(A, B);

    case 13
      A = read_matrix('A');
      B = read_matrix('B');

    case 0
      running = false;
      disp('Goodbye!');

    otherwise
      disp('Invalid option. Please choose a number from the menu.');
  end
end
