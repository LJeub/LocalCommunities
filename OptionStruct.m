classdef OptionStruct < matlab.mixin.Copyable
    % Class for option parsing based on structs.
    % Construct with OptionStruct(options), where options is either:
    %       a list of possible options
    %       key value pairs (listing option names and their defaults
    %       a struct
    %       an OptionStruct (copies the option struct)
    % new options can only be added at construction time
    %
    % methods: set: set option values (does not create new fields and warns
    % about non-existing options)
    %
    % properties: options: used to get a list of options
    %
    % options are referred to using struct syntax, i.e. OptionStruct.option

% Version:
% Date:
% Author:
% Email:
    
    properties (Hidden)
        option_struct=struct([]);
    end
    
    properties (Dependent=true, SetAccess=private)
        options
    end
    
    methods
        %constructor (pass struct with defaults or list of options)
        function obj=OptionStruct(varargin)
            if nargin>0
                input=varargin;
                %unpack nested cells (usually created by passing varargins)
                while iscell(input)&&length(input)==1
                    input=input{1};
                end
                if length(input)==1
                    if isstruct(input)
                        obj.options=fieldnames(varargin{1});
                        obj.set(varargin{1});
                    elseif isa(input,'OptionStruct')
                        obj=copy(input);
                    else
                        obj.options=input;
                    end
                else
                    if iscellstr(input)
                        obj.options=input;
                    else
                        if ~mod(length(input),2)
                            for i=1:2:length(input)
                                if isvarname(lower(input{i}))
                                    obj.options=lower(input{i});
                                    obj.option_struct.(lower(varargin{i}))=varargin{i+1};
                                else
                                    error('%s is not a valid option name',input{i})
                                end
                            end
                        else
                            error('option list has odd length')
                        end
                    end
                end
            end
        end
        
        
        %subscripted assignment
        function obj=subsasgn(obj, S, value)
            if isequal(S(1).type,'.')
                if ismember(S(1).subs,properties(obj))||ismember(S(1).subs,methods(obj))
                    obj=builtin('subsasgn',obj,S,value);
                    return;
                end
                
                S(1).subs=lower(S(1).subs);
                if obj.isfield(S(1).subs)
                    obj.option_struct=builtin('subsasgn',obj.option_struct,S,value);
                else
                    error('option %s does not exist',S.subs);
                end
            else
                error('subscripted assignment %s undefined',S.type);
            end
        end
        
        %subscripted reference
        function val=subsref(obj,S)
            if isequal(S(1).type,'.')
                try
                    val=builtin('subsref',obj,S);
                    return;
                catch err
                    if strcmp(err.identifier,'MATLAB:noSuchMethodOrField')
                        S(1).subs=lower(S(1).subs);
                        if obj.isfield(S(1).subs)
                            val=builtin('subsref',obj.option_struct,S);
                        else
                            error('option %s does not exist',S.subs);
                        end
                    else
                        rethrow(err);
                    end
                end
            else
                error('subscripted reference %s is undefined',S(1).type);
            end
        end
        
        %nice display of options
        function disp(obj)
            disp('option struct with fields:')
            disp(obj.option_struct)
        end
        
        %set and get possible option names
        function set.options(obj,fieldnames)
            if ischar(fieldnames)
                fieldnames={fieldnames};
            end
            if ~iscellstr(fieldnames)
                error('cannot parse fieldnames');
            end
            fieldnames=lower(fieldnames);
            for i=1:length(fieldnames)
                obj.option_struct(1).(fieldnames{i})=[];
            end
        end
        function opt=get.options(obj)
            opt=fieldnames(obj.option_struct);
        end
        
        %set options (takes a 'struct' or 'key','value' pairs)
        function obj=set(obj,varargin)
            input=varargin;
            while iscell(input)&&length(input)==1
                input=input{1};
            end
            if length(input)==1
                if isstruct(input)
                    fields=fieldnames(input);
                    for i=1:length(fields)
                        if obj.isfield(fields{i})
                            obj.option_struct.(lower(fields{i}))=input.(fields{i});
                        else
                            error('option ''%s'' does not exist',fields{i});
                        end
                    end
                else
                    error('need input struct');
                end
            else
                obj.parse_option_list(input);
            end
        end
        
        function is_opt=isfield(obj,fieldname)
            is_opt=isfield(obj.option_struct,lower(fieldname));
        end
        
        function is_set=isset(obj,fieldname)
            fieldname=lower(fieldname);
            if ~iscell(fieldname)
                fieldname={fieldname};
            end
            opt=obj.isfield(fieldname);
            if any(~opt)
                error('option ''%s'' does not exist \n',fieldname{~opt});
            end
            
            is_set=false(size(fieldname));
            for i=1:length(fieldname)
                is_set(i)=~isempty(obj.option_struct.(fieldname{i}));
            end
        end
        
    end
    
    methods (Access=private)
        function parse_option_list(obj,list)
            if iscell(list)
                if ~mod(length(list),2)
                    for i=1:2:length(list)
                        if obj.isfield(list{i})
                            obj.option_struct.(lower(list{i}))=list{i+1};
                        else
                            optstr=sprintf('%s \n',obj.options{:});
                            error('option ''%s'' does not exist, valid options are %s',list{i},optstr);
                        end
                    end
                else
                    error('option list has odd length')
                end
            else
                error('option list has to be a cell array')
            end
        end
    end
    
end

