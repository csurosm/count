package ca.umontreal.iro.banality;

import java.util.ArrayList;

/**
 *
 * Simple parsing for command-line switches
 *
 * @author csuros
 */
public class CommandLineSwitches
{
    private ArrayList<String> arguments;

    public CommandLineSwitches()
    {
        arguments = new ArrayList<String>();
    }

    public CommandLineSwitches(String[] args)
    {
        arguments = new ArrayList<String>();
        for (String a: args)
            arguments.add(a);
    }

    public void clear()
    {
        arguments.clear();
    }

    public String nextSwitch()
    {
        if (arguments.size()>0 && arguments.get(0).startsWith("-"))
            return arguments.get(0);
        else
            return null;
    }

    public String[] getArguments()
    {
        return arguments.toArray(new String[0]);
    }

    public <U> void inspect(U u)
    {
        System.out.println("U: " + u.getClass().getName());
    }


    public <T> T parse(String switch_name, TypeParser<? extends T> parser)
    {
        int idx = arguments.indexOf(switch_name);
        if (idx==-1)
            return parser.parse(null);
        {
            arguments.remove(idx);
            if (arguments.size()==idx)
                throw new IllegalArgumentException("Switch "+switch_name+" needs a value ");
            return parser.parse(arguments.remove(idx));
        }
    }

    public Object parse(String switch_name, SwitchParser parser)
    {
        int idx = arguments.indexOf(switch_name);
        if (idx==-1)
            return parser.parse(null);
        else
        {
            arguments.remove(idx);
            if (arguments.size()==idx)
                throw new IllegalArgumentException("Switch "+switch_name+" needs a value ");
            return parser.parse(arguments.remove(idx));
        }
    }

    public Object parse(String switch_name, Object default_value)
    {
        return parse(switch_name, new ObjectParser(default_value));
    }

    public boolean hasSwitch(String switch_name)
    {
        int idx = arguments.indexOf(switch_name);
        if (idx==-1)
            return false;
        else
        {
            arguments.remove(idx);
            return true;
        }
    }

    public boolean parseBoolean(String switch_name, boolean default_value)
    {
        return (Boolean) parse(switch_name, new BooleanParser(default_value));
    }

    public double parseDouble(String switch_name, double default_value)
    {
        return (Double) parse(switch_name, new DoubleParser(default_value));
    }

    public String parseString(String switch_name, String default_value)
    {
        return (String) parse(switch_name, new StringParser(default_value));
    }

    public int parseInt(String switch_name, int default_value)
    {
        return (Integer) parse(switch_name, new IntegerParser(default_value));
    }

    public long parseLong(String switch_name, long default_value)
    {
        return (Long) parse(switch_name, new LongParser(default_value));
    }


    private static final SwitchParser IDENTITY_PARSER = new StringParser(null);

    public static class ObjectParser extends TypeParser<Object>
    {
        public ObjectParser(Object o)
        {
            super(o);
        }

        @Override
        public Object convert(String val)
        {
            return val;
        }
    }

    public static class IntegerParser extends TypeParser<Integer>
    {
        public IntegerParser(Integer i)
        {
            super(i);
        }

        @Override
        public Integer convert(String val)
        {
            return Integer.parseInt(val);
        }
    }

    public static class LongParser extends TypeParser<Long>
    {
        public LongParser(Long l)
        {
            super(l);
        }

        @Override
        public Long convert(String val)
        {
            return Long.parseLong(val);
        }
    }

    public static class StringParser extends TypeParser<String>
    {
        public StringParser(String s)
        {
            super(s);
        }

        @Override
        protected String convert(String value)
        {
            return value;
        }
    }

    public static class BooleanParser extends TypeParser<Boolean>
    {
        public BooleanParser(Boolean default_value)
        {
            super(default_value);
        }

        @Override
        protected Boolean convert(String value)
        {
            return "true".equals(value);
        }
    }

    public static class DoubleParser extends TypeParser<Double>
    {
        public DoubleParser(Double default_value)
        {
            super(default_value);
        }

        @Override
        protected Double convert(String value)
        {
            return Double.parseDouble(value);
        }
    }

    
    public static abstract class TypeParser<T> implements SwitchParser
    {
        private T default_value;
        public TypeParser(T default_value)
        {
            this.default_value = default_value;
        }
        @Override
        public T parse(String switch_value)
        {
            if (switch_value == null)
                return default_value;
            else
                return convert(switch_value);
        }
        
        protected abstract T convert(String value);
    }

    public interface SwitchParser<T>
    {
        public T parse(String switch_value);
    }
}
