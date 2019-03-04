% colours pallete definition based on vibrant Paul Tol colors

switch pallete
    case 'bright'
        % RGB Bright qualitative colour scheme
        BLUE = [68,119,170];
        CYAN = [102,204,238];
        GREEN = [34,136,51];
        YELLOW = [204,187,68];
        RED = [238,102,119];
        PURPLE = [170,51,119];
        GREY = [187,187,187];
        COLOURS = {BLUE,RED,GREEN,YELLOW,CYAN,PURPLE,GREY};
    case 'high-contrst'
        % RGB High-contrast qualitative colour scheme (good greyscale reproduction)
        WHITE = [255,255,255];
        YELLOW = [221,170,51];
        RED = [187,85,102];
        BLUE = [0,68,136];
        BLACK = [0,0,0];
        COLOURS ={BLUE,YELLOW,RED};
    case 'vibrant'
        % RGB Vibrant qualitative colour scheme
        BLUE = [0,119,187];
        CYAN = [51,187,238];
        TEAL = [0,153,136];
        ORANGE = [238,119,51];
        RED = [204,51,17];
        MAGENTA = [238,51,119];
        GREY = [187,187,187];
        COLOURS = {ORANGE,BLUE,CYAN,MAGENTA,RED,TEAL,GREY};
    case 'muted'
        INDIGO = [51,34,136];
        CYAN = [136,204,238];
        TEAL = [68,170,153];
        GREEN = [17,119,51];
        OLIVE = [153,153,51];
        SAND = [221,204,119];
        ROSE = [204,102,119];
        WINE = [136,34,85];
        PURPLE = [170,68,153];
        PALE_GREY = [221,221,221];
        COLOURS = {ROSE,INDIGO,SAND,GREEN,CYAN,WINE,TEAL,OLIVE,PURPLE};
    case 'light'
        % RGB Light qualitative colour scheme
        LIGHT_BLUE = [119,170,221];
        LIGHT_CYAN = [153,221,255];
        MINT = [68,187,153];
        PEAR = [187,204,51];
        OLIVE = [170,170,0];
        LIGHT_YELLOW = [238,221,136];
        ORANGE = [238,136,102];
        PINK = [255,170,187];
        PALE_GREY = [221,221,221];
        COLOURS = {LIGHT_BLUE,ORANGE,LIGHT_YELLOW,PINK,LIGHT_CYAN,MINT,PEAR,OLIVE,PALE_GREY};
end

% Diverging colour schemes
% Diverging schemes are for ordered data between two extremes 
% where the midpoint is important
% colormaps that also works in colour-blind vision
map_sunset = [54,75,154;
               74,123,183;
               110,166,205;
               152,202,225;
               194,228,239;
               234,236,204;
               254,218,139;
               253,179,102;
               246,126,75;
               221,61,45;
               165,0,38];
map_burd = [33,102,172;
             67,147,195;
             146,197,222;
             209,229,240;
             247,247,247;
             253,219,199;
             244,165,130;
             214,96,77;
             178,24,43];

N=256; % number of rows in colormap
x_interp = 0:1:N-1; %
x1 = linspace(0,N-1,length(map_sunset));
map_sunset_interp = interp1(x1,map_sunset,x_interp','linear','extrap');
x2 = linspace(0,N-1,length(map_burd));
map_burd_interp = interp1(x2,map_burd,x_interp','linear','extrap');
% Sequential colour schemes
% Sequential schemes are for ordered data from low to high
map_brewer = [255,255,229;
             255,247,188;
             254,227,145;
             254,149,79;
             251,154,41;
             236,112,20;
             204,76,2;
             153,52,4;
             102,37,6];
map_iridescent = [254,251,233;
                  254,247,213;
                  245,243,193;
                  234,240,181;
                  221,236,191;
                  208,231,202;
                  194,227,210;
                  181,221,216;
                  168,216,220;
                  155,210,225;
                  129,196,231;
                  123,188,231;
                  126,178,228;
                  136,165,221;
                  147,152,210;
                  155,138,196;
                  157,125,178;
                  154,112,158;
                  144,99,136;
                  128,87,112;
                  104,73,87;
                  70,53,58];
x3 = linspace(0,N-1,length(map_brewer));
map_brewer_interp = interp1(x3,map_brewer,x_interp','linear','extrap');
x4 = linspace(0,N-1,length(map_iridescent));
map_iridescent_interp = interp1(x4,map_iridescent,x_interp','linear','extrap');   

map_discrete_rainbow = [209,187,215;
                        174,118,163;
                        136,46,114;
                        25,101,176;
                        82,137,199;
                        123,175,222;
                        78,178,101;
                        144,201,135;
                        202,224,171;
                        247,240,86;
                        246,193,65;
                        241,147,45;
                        232,96,28;
                        220,5,12];