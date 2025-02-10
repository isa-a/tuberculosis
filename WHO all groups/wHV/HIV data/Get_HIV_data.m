clear all; 

countries = {'Congo','Ethiopia','Kenya','Malawi','Mozambique','Nigeria','South Africa','United Republic of Tanzania','Uganda','Zambia','Zimbabwe'};

% --- Get estimates on new infections -------------------------------------

C    = readcell('New_infections.xlsx');
rows = [];
for ir = 1:length(countries)
   rows(ir) = find(strcmp(C(:,1),countries{ir})); 
end
dat  = C(rows,2:end);

mat = [];
for ic = 1:size(dat,2)
    for ir = 1:size(dat,1)
        txt  = dat{ir,ic};
        txt2 = txt(find(~isspace(txt)));
        txt3 = strrep(txt2,'[','-');
        txt4 = txt3(1:end-1);
        txt5 = strsplit(txt4,'-');
        mat(ir,ic,:) = [str2num(txt5{2}), str2num(txt5{1}), str2num(txt5{3})];
    end
end

mat = squeeze(sum(mat,1));

% Assume HIV burden scaled up linearly from 1980 to first data point in
% 1990
yr0  = 1980;
Ys   = [zeros(1,3); mat(1,:)];
Xs   = [1980 1990];
ms   = diff(Ys,1)/diff(Xs);
cs   = Ys(1)-ms*Xs(1);
xpts = [Xs(1):1:Xs(2)]';
ypts = xpts.*ms + cs;

HIV_incd = [ypts; mat];


% --- Get estimates on ART coverage ---------------------------------------

C    = readcell('ART coverage.xlsx');
rows = [];
for ir = 1:length(countries)
   rows(ir) = find(strcmp(C(:,1),countries{ir})); 
end
dat  = C(rows,2:end);

mat = [];
for ic = 1:size(dat,2)
    for ir = 1:size(dat,1)
        txt  = dat{ir,ic};
        txt2 = txt(find(~isspace(txt)));
        txt3 = strrep(txt2,'[','-');
        txt4 = txt3(1:end-1);
        txt5 = strsplit(txt4,'-');
        mat(ir,ic,:) = [str2num(txt5{2}), str2num(txt5{1}), str2num(txt5{3})];
    end
end
% figure; plot(mat);
% yl = ylim; yl(1) = 0; ylim(yl);


ARTcovg_2019 = mat(end,:);

% Do linear extrapolation backwards to figure out when ART started
Ys     = mat([1,end],2);
Xs     = [2010 2019];
ms     = diff(Ys)/diff(Xs);
cs     = Ys(1)-ms*Xs(1);
ART_start = floor(-cs/ms);


% --- Get estimates on HIV prevalence -------------------------------------

C = readcell('HIV prevalence.xlsx');
row = find(strcmp(C(:,1),country));
vec = C(row,2:end);

mat = [];
for iv = 1:length(vec)
    txt  = vec{iv};
    txt2 = txt(find(~isspace(txt)));
    txt3 = strrep(txt2,'[','-');
    txt4 = txt3(1:end-1);
    txt5 = strsplit(txt4,'-');
    mat(iv,:) = [str2num(txt5{2}), str2num(txt5{1}), str2num(txt5{3})];
end
% figure; plot(mat);

% Get the last available year of data
HIVprev_2019 = mat(end,:);


% --- Save all under the relevant filename --------------------------------

fname = ['HIV_inputs_',country(find(~isspace(country)))];
save(fname,'HIV_incd','ART_start','ARTcovg_2019','HIVprev_2019');
