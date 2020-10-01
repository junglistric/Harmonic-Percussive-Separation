%harmonic suppression processing
%function output = HPS(filename,win_size,hop_size,medfilt_len,outfile)
function output1 = HPS(filename,win_size,hop_size,medfilt_len,outfile)

%read signals
[wave1,fs]=audioread(filename);

%starting sample
startingSamp = 1;

%create the sine window function
win = sin(pi/win_size*(1:win_size)');

%create a half cosine window
coswin = win(win_size/2+1:end);

%create a half sine window
sinwin = win(1:win_size/2);

%number of loops + 1
numblocks = floor(length(wave1)/hop_size);

%determine how many samples to pad
pad = length(wave1)-(numblocks*hop_size);

%pad the end of the signals with zeros
wave1 = [wave1;zeros(pad,2)];

%output vector
output1 = zeros(size(wave1));

%first block
%create a vector containing a block of signal
chunka = wave1(startingSamp : startingSamp + hop_size - 1, 1);
    
%apply window
winchunka = coswin .* chunka;
    
%take fft
WINCHUNKA = fft(winchunka);

%cartesian to polar
[th,r] = cart2pol(real(WINCHUNKA),imag(WINCHUNKA));
r_out = zeros(length(r),1);

%freq domain median filter lengths
if mod(medfilt_len,2) == 1
    
    L1 = ((medfilt_len-1)/2);
    L2 = ((medfilt_len-1)/2);
    
elseif mod(medfilt_len,2) == 0
   
    L1 = ceil((medfilt_len-1)/2);
    L2 = floor((medfilt_len-1)/2);
    
end

%median filter
for m = 1:L1

    r_out(m) = median(r(1:m+L2));

end

for m = L1+1:length(r)/2-L2

    r_out(m) = median(r(m-L1:m+L2));

end

for m = length(r)/2-L2+1:length(r)/2

    r_out(m) = median(r(m-L1:length(r)/2));

end

%mirror until length of r
R = [r_out(1:length(r_out)/2);r_out(length(r_out)/2:-1:1)];

%back to cartesian
[re, im] = pol2cart(th,R);
WINCHUNKA = re + i*im;

%take inverse fft
winchunka = real(ifft(WINCHUNKA));
    
%apply window
winchunka = coswin .* winchunka;
    
%overlap and add
output1(startingSamp:startingSamp + hop_size - 1, 1) = output1(startingSamp:startingSamp + hop_size - 1, 1) + winchunka;

%iterate over blocks
for j = 1:numblocks - 1

    %create a vector containing a block of signal
    chunk1a = wave1(startingSamp : startingSamp + hop_size - 1, 1);
    
    %retain the starting sample position
    startingSamp1 = startingSamp;
        
    %set the next starting sample position
    startingSamp = (startingSamp + win_size) - hop_size;
    
    %segment second block
    chunk1b = wave1(startingSamp : startingSamp + hop_size - 1, 1);

    %concatenate chunks
    chunk1 = [chunk1a;chunk1b];
    
    %apply window
    winchunk1 = win .* chunk1;

    %take fft
    WINCHUNK1 = fft(winchunk1);

    %cartesian to polar
    [th,r] = cart2pol(real(WINCHUNK1),imag(WINCHUNK1));
    r_out = zeros(length(r),1);

    %median filter
    for m = 1:L1

        r_out(m) = median(r(1:m+L2));

    end

    for m = L1+1:length(r)/2-L2

        r_out(m) = median(r(m-L1:m+L2));

    end

    for m = length(r)/2-L2+1:length(r)/2

        r_out(m) = median(r(m-L1:length(r)/2));

    end

    %mirror until length of r
    R = [r_out(1:length(r_out)/2);r_out(length(r_out)/2:-1:1)];

    %back to cartesian
    [re, im] = pol2cart(th,R);
    WINCHUNKA = re + i*im;

    %take inverse fft
    winchunk1 = real(ifft(WINCHUNKA));
    
    %apply window
    winchunk1 = win .* winchunk1;
    
    %overlap and add
    output1(startingSamp1:startingSamp + hop_size - 1, 1) = output1(startingSamp1:startingSamp + hop_size - 1, 1) + winchunk1;

end

%final block
%create a vector containing a block of audio
chunka = wave1(startingSamp : startingSamp + hop_size - 1, 1);
    
%retain the starting sample position
startingSamp1 = startingSamp;
        
%set the next starting sample position
startingSamp = (startingSamp + win_size) - hop_size;
    
%apply window
winchunka = sinwin .* chunka;
    
%take fft
WINCHUNKA = fft(winchunka);
   
%cartesian to polar
[th,r] = cart2pol(real(WINCHUNKA),imag(WINCHUNKA));
r_out = zeros(length(r),1);

%median filter
    for m = 1:L1

        r_out(m) = median(r(1:m+L2));

    end

    for m = L1+1:length(r)/2-L2

        r_out(m) = median(r(m-L1:m+L2));

    end

    for m = length(r)/2-L2+1:length(r)/2

        r_out(m) = median(r(m-L1:length(r)/2));

    end

    %mirror until length of r
    R = [r_out(1:length(r_out)/2);r_out(length(r_out)/2:-1:1)];

    %back to cartesian
    [re, im] = pol2cart(th,R);
    WINCHUNKA = re + i*im;

%take inverse fft
winchunka = real(ifft(WINCHUNKA));
    
%apply window
winchunka = sinwin .* winchunka;
    
%overlap and add
output1(startingSamp1:startingSamp - 1, 1) = output1(startingSamp1:startingSamp - 1, 1) + winchunka;

%trim outputs
output1=output1(1:end-pad, :);

%output the result
audiowrite(outfile,output1,fs);

end