classdef PlotStyle
    %PLOTSTYLE  IEEE 论文级配色/字体/样式工具包
    %   全局配色方案 + 一键 axes 样式设置, 适用于所有仿真绘图脚本.
    %
    %   用法:
    %     ps = PlotStyle;                       % 实例化
    %     plot(x, y, '-o', 'Color', ps.Blue);   % 取色
    %     ps.apply(gca);                        % 一键设置当前 axes
    %     ps.apply(gca, 10);                    % 指定字号
    %     fig = ps.newFigure(720, 500);         % 创建标准尺寸白底图窗
    %
    %   配色一览:
    %     Blue / Red / Green / Purple / Orange / Gray / Black
    %     BlueDark / RedDark / GreenDark / PurpleDark  (深色变体)
    %     BlueLight / RedLight / GreenLight             (浅色填充)

    %% ================ 主色板 ================
    properties (Constant)
        % --- 主色 (实线/Marker) ---
        Blue    = [0.00 0.45 0.74]
        Red     = [0.85 0.33 0.10]
        Green   = [0.47 0.67 0.19]
        Purple  = [0.49 0.18 0.56]
        Orange  = [0.93 0.69 0.13]
        Gray    = [0.50 0.50 0.50]
        Black   = [0.00 0.00 0.00]

        % --- 深色变体 (图例/标注) ---
        BlueDark   = [0.00 0.30 0.55]
        RedDark    = [0.65 0.20 0.05]
        GreenDark  = [0.30 0.50 0.10]
        PurpleDark = [0.35 0.10 0.42]

        % --- 浅色变体 (填充/阴影/置信区间) ---
        BlueLight  = [0.68 0.82 0.93]
        RedLight   = [0.96 0.76 0.65]
        GreenLight = [0.80 0.90 0.65]

        % --- 字体 ---
        FontName = 'Times New Roman'
        FontSize = 10              % 默认 axes 字号
        TitleSize = 12             % sgtitle 字号
    end

    %% ================ 有序取色 ================
    methods
        function c = color(obj, idx)
            %COLOR  按序号取色, 循环使用
            %   c = ps.color(1)  → Blue
            %   c = ps.color(5)  → Orange
            palette = [obj.Blue; obj.Red; obj.Green; ...
                       obj.Purple; obj.Orange; obj.Gray; obj.Black];
            c = palette(mod(idx - 1, size(palette, 1)) + 1, :);
        end
    end

    %% ================ 样式设置 ================
    methods
        function apply(obj, ax, fontSize)
            %APPLY  一键设置 axes 为 IEEE 论文样式
            %   ps.apply(gca)       使用默认字号
            %   ps.apply(gca, 11)   指定字号
            if nargin < 3
                fontSize = obj.FontSize;
            end

            set(ax, ...
                'FontName',     obj.FontName, ...
                'FontSize',     fontSize, ...
                'TickDir',      'in', ...
                'TickLength',   [0.015 0.015], ...
                'Box',          'on', ...
                'LineWidth',    0.8, ...
                'XMinorTick',   'on', ...
                'YMinorTick',   'on', ...
                'TickLabelInterpreter', 'latex');

            grid(ax, 'on');
            ax.GridAlpha      = 0.3;
            ax.MinorGridAlpha = 0.15;
        end

        function fig = newFigure(obj, w, h)
            %NEWFIGURE  创建白底图窗, 默认 560×420 (IEEE 半栏)
            %   fig = ps.newFigure;
            %   fig = ps.newFigure(800, 600);
            if nargin < 2, w = 560; end
            if nargin < 3, h = 420; end
            fig = figure('Position', [100 100 w h], 'Color', 'w');
        end

        function applyLegend(obj, lg, fontSize)
            %APPLYLEGEND  设置图例样式
            %   lg = legend(...); ps.applyLegend(lg);
            if nargin < 3
                fontSize = obj.FontSize - 1;
            end
            set(lg, ...
                'FontName',    obj.FontName, ...
                'FontSize',    fontSize, ...
                'Interpreter', 'latex', ...
                'Box',         'off');
        end
    end

    %% ================ 快捷绘图 ================
    methods
        function h = plotLine(obj, ax, x, y, colorIdx, varargin)
            %PLOTLINE  带自动配色的实线绘制
            %   h = ps.plotLine(gca, x, y, 1);             → 蓝色实线
            %   h = ps.plotLine(gca, x, y, 2, '--', 'o');  → 红色虚线圆圈
            %
            %   可选参数: lineStyle, marker (按位置传入)
            lineStyle = '-';
            marker    = 'none';
            if numel(varargin) >= 1, lineStyle = varargin{1}; end
            if numel(varargin) >= 2, marker    = varargin{2}; end

            c = obj.color(colorIdx);
            h = plot(ax, x, y, ...
                'LineStyle',       lineStyle, ...
                'Marker',          marker, ...
                'Color',           c, ...
                'MarkerFaceColor', c, ...
                'MarkerSize',      5, ...
                'LineWidth',       1.2);
        end
    end

    %% ================ 静态工具 ================
    methods (Static)
        function demo()
            %DEMO  展示完整色板
            ps = PlotStyle;
            fig = ps.newFigure(600, 300);

            names = {'Blue','Red','Green','Purple','Orange','Gray','Black'};
            for ii = 1:7
                x = linspace(0, 2*pi, 100);
                y = sin(x + ii*0.5) + ii*1.5;
                ps.plotLine(gca, x, y, ii);
                hold on;
            end
            hold off;

            lg = legend(names, 'Location', 'eastoutside');
            ps.applyLegend(lg);
            ps.apply(gca);
            title('PlotStyle Color Palette Demo', ...
                'Interpreter', 'latex', 'FontSize', ps.TitleSize);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel('$y$', 'Interpreter', 'latex');

            fprintf('  PlotStyle demo 完成.\n');
        end
    end
end