package count.machine;

import java.util.Map;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;

/**
 * Common foundation class for Machines. 
 * 
 * @author csuros
 *
 */
public abstract class Machine 
{
	static final boolean PRINT_LIKELIHOODS = true;
	
	Machine()
	{
		Class<? extends Machine> whatIam = this.getClass();
		Type t = typeOf(whatIam);
		long count = machine_counts.get(t);
		++ count;
		this.machine_id = t.toString().substring(0,1)+Long.toString(count);
		machine_counts.put(t,  count);
	}
	
	final String machine_id;
	
	static enum Type {Other, Node, Insert, Mutate};
	final static Map<Type, Long> machine_counts;

	static {
		machine_counts = new EnumMap<>(Type.class);
		for (Type t: Type.values())
		{
			machine_counts.put(t,  0L);
		}
	}
	static Type typeOf(Class<? extends Machine> C)
	{
		
		
		if (C==Node.class || C==Differ.class)
		{
			return Type.Node;
		} else if (C==Insert.class)
		{
			return Type.Insert;
		} else if (C == Mutate.class || C==Root.class)
		{
			return Type.Mutate;
		}
		return Type.Other;
	}
	
	
	abstract void computeReadLikelihoods();
	abstract void computeWriteProbabilities();
	
	/**
	 * Whether this is a terminal machine that does not write
	 * @return true for leaf nodes containing the observed sequences
	 */
	abstract boolean isLeaf();
	/**
	 * Whether this is a root machine with no upstream connection
	 * @return true for the root's insert machine
	 */
	abstract boolean isRoot();
	/**
	 * Postorder traversal of the connected downstream machines. The default implementation
	 * for a leaf machine with no downstream connections adds only this to the machine list.    
	 * 
	 * @param <T> class for the machines to be included
	 * @param machines null for first call
	 * @param machineClass (Machine or any of its subclasses)
	 * @return a list of machines with leaves first and root last. 
	 */
	<T extends Machine> List<T> postOrder(List<T> machines, Class<T> machineClass)
	{
		if (machines==null)
			machines = new ArrayList<>();

		if (machineClass.isInstance(this))
			machines.add((T) this); // meaningless type cast (T), but we are sure of the instance
		return machines;
	}
	
	
	
	@Override
	public String toString()
	{
		return machine_id;
	}

}
