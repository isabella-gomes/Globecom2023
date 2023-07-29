% Generate uncorrelated communication symbols for MIMO using spatial multiplexing

function [transmitSymbols,waveforms]=signalgeneration(nTx,num_waveforms)
% Set the number of transmit and receive antennas
%nTx = 4;
%nRx = 4;

% Set the number of subcarriers and the modulation order
nSubcarriers = 128;
modulationOrder = 16;

% Generate random data symbols
dataSymbols = randi([0, modulationOrder-1], nTx, nSubcarriers);

% Assign phase and amplitude values to each data symbol
phase = exp(1i*2*pi*rand(nTx, nSubcarriers));
amplitude = sqrt(1/nTx)*ones(nTx, nSubcarriers);

% Modulate each data symbol using QAM
modulatedSymbols = qammod(dataSymbols, modulationOrder, 'UnitAveragePower', true);

% Apply phase and amplitude values to the modulated symbols
transmitSymbols = bsxfun(@times, modulatedSymbols, amplitude) .* phase;


% Generate individual radar waveforms with pseudo random coding

% Define number of waveforms and waveform length
%num_waveforms = 4;
waveform_length = 128;

% Generate pseudo-random codes for each waveform
codes = zeros(num_waveforms, waveform_length);
for i = 1:num_waveforms
    codes(i,:) = sign(randn(1, waveform_length));
end

% Calculate cross-correlation matrix between codes
%corr_matrix = codes * codes';

% Generate uncorrelated waveforms using Hadamard matrix
hadamard_matrix = hadamard(waveform_length);
waveforms = zeros(num_waveforms, waveform_length);
for i = 1:num_waveforms
    waveforms(i,:) = hadamard_matrix(i,:) * sqrt(waveform_length);
end

% Multiply each waveform with its corresponding code and normalize
for i = 1:num_waveforms
    waveforms(i,:) = waveforms(i,:) .* codes(i,:) / norm(codes(i,:));
end
end