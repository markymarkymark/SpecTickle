function [fid,hdr,t,hz,ppm] = siemens_dicom_mrs(hdr,opt)
% spec_xaxis - calculate time and frequency x-axis vectors for spectroscopy
% 
% Syntax:
%   [fid,hdr,t,hz,ppm] = read_siemens_mrs(hdr,opt)
%
% Created: Mark A. Elliott, PhD
%   melliott@upenn.edu
%   https://www.med.upenn.edu/CAMIPM/mark-elliott.html

fid  = [];

% --- Prompt user for a file ---
if (nargin < 1)
    hdr = dicom_header();  
    if (isempty(hdr)); return; end
end

% --- Check options ---
if (nargin < 2), opt.dummy = 1; end
opt.verbose      = checkstruct(opt,'verbose',1);    
opt.read_skippts = checkstruct(opt,'read_skippts',0);   % discard initial bad data points 
if (opt.verbose)
    fprintf(1,'Reading Siemens Dicom MRS file: %s\n',hdr.Filename);
    if (opt.read_skippts > 0), fprintf(1,'    Skipping %1d initial Dicom FID points\n',opt.read_skippts); end
end

% --- Read spectro data from Siemens header ---
if (~isfield(hdr,'Private_7fe1_1010')), fprintf(2,'ERROR: File is not a Siemens MRS Dicom\n'); return; end
nums = hdr.Private_7fe1_1010;
data = typecast(nums,'single');
N    = size(data,1)/2;
even = 1:2:2*N;
fid  = complex(data(even),data(even+1));

if (opt.read_skippts > 0) % NOTE: this isn't coded for CSI data!!
    fid = circshift(fid,-opt.read_skippts);
    fid(end-opt.read_skippts+1:end) = fid(end-opt.read_skippts);
    fid = fid * exp(complex(0,-angle(fid(1)))); % adjust constant phase
end

% --- data info - need to get from Siemens Private header ---
[~,f0,hdr] = dicom_get_header(hdr,'ImagingFrequency');
[~,~,hdr] = dicom_get_header(hdr,'PixelBandwidth');
[~,~,hdr] = dicom_get_header(hdr,'RepetitionTime');
[~,~,hdr] = dicom_get_header(hdr,'EchoTime');
[~,dwell,hdr] = dicom_get_header(hdr,'RealDwellTime'); % BW is 1/(hdr.dwell * 1e-9) meaining dwell is in nsec!
%[~,shadow,hdr] = dicom_get_header(hdr,'MrPhoenixProtocol'); % look for string 'sSpecPara.ucRemoveOversampling' in shadow 

% --- Add params to hdr ---
hdr.MRS_f0 = f0;
hdr.MRS_N  = N;
hdr.MRS_BW = 1/(dwell * 1e-9);	

% --- return additional items ---
[t,hz,ppm] = spec_xaxis(N,hdr.MRS_BW,f0);

% --- CSI data ---
if (hdr.is_csi)
    [~,nx] = dicom_get_header(hdr,'SpectroscopyAcquisitionPhaseRows');
    [~,ny] = dicom_get_header(hdr,'SpectroscopyAcquisitionPhaseColumns');
%    [~,nz] = dicom_get_header(hdr,'SpectroscopyAcquisitionOut_of_planePhaseSteps');
    [~,nz] = dicom_get_header(hdr,'NumberOfFrames');
%    [~,nz] = dicom_get_header(hdr,'SpectroscopyAcquisitionDataColumns');
    [~,nt] = dicom_get_header(hdr,'DataPointColumns');
    
    hdr.CSI_Nx = nx;
    hdr.CSI_Ny = ny;
    hdr.CSI_Nz = nz;
    hdr.CSI_Nt = nt;
    hdr.CSI_BW = 1/(dwell * 1e-9);
    fid = reshape(fid,nt,nx,ny,nz);
end
return;
