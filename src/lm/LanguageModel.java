package lm;

import java.util.List;

/**
 * @author jda
 */
public interface LanguageModel<T> {

    public T sample(List<T> context);

    public double score(T entry, List<T> context);

}
