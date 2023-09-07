function [s, w] = runFslCommand(rap,fslcmd,ENV,varargin)

% Setup
fslsetup = deblank(rap.directoryconventions.fslsetup);
if not(isempty(fslsetup))
    if ~startsWith(fslsetup,'. '), fslsetup = ['. ' fslsetup]; end
    if ~endsWith(fslsetup,';'), fslsetup=[fslsetup ';']; end
end

% Parse (e.g., when called from runPyCommand)
indShPfx  = find(strcmp(varargin,'shellprefix'));
if ~isempty(indShPfx)
    fslsetup = [varargin{indShPfx+1} fslsetup];
    varargin(indShPfx:indShPfx+1) = [];
end

if nargin < 3, ENV = {}; end
if nargin < 5, varargin = {}; end

global reproacache
if ~reproacache.isKey('toolbox.spm'), logging.error('SPM is not found'); end
SPMtool = reproacache('toolbox.spm');

ENV = vertcat(ENV,{...
    'FSLOUTPUTTYPE', rap.directoryconventions.fsloutputtype; ...
    'MATLABPATH', SPMtool.toolPath;...
    });

[s, w] = shell(fslcmd,'shellprefix',fslsetup,'environment',ENV,varargin{:});
